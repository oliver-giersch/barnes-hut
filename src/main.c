// Required for `pthread_barrier`.
#ifdef __linux
#define _XOPEN_SOURCE 600
#endif // __linux

#include <stdatomic.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <errno.h>
#include <pthread.h>
#include <string.h>
#include <time.h>

#include <sys/mman.h>

#include "barnes-hut/arena.h"
#include "barnes-hut/common.h"
#include "barnes-hut/options.h"
#include "barnes-hut/phys.h"
#ifdef RENDER
#include "barnes-hut/render.h"
#endif // RENDER

// The per-thread simulation state.
struct thread_state {
	// The thread's ID.
	unsigned id;
	// The thread's local copy of the global particle list.
	struct accel_particle *particles;
	// The thread's assigned particle slice.
	struct particle_slice slice;
	// The particle space radius to use for the next simulation step.
	//
	// After completion of a step, the thread sets this field to the greatest
	// distance from center of any calculated particle position in the latest
	// iteration.
	//
	// The main thread then aggregates all thread values of this field and
	// calculates (maximum) the radius for the next iteration.
	//
	// Access to this field must be synchronized using `barrier`.
	float radius;
} aligned(64);

// The global thread error flag.
static atomic_int thread_errno = 0;
// The global thread synchronization barrier.
static pthread_barrier_t barrier;
// The globally shared and synchronized region of all simulated particles.
static struct accel_particle *particles;
// The memory arena for octant allocation.
static struct arena arena;
// The globally shared and synchronized tree of particles.
//
// Access to the tree must be synchronized using `barrier`.
static struct particle_tree tree;
// The TLS holding the state of all threads.
static struct threads {
	unsigned len;
	struct thread_state states[];
} *tls = NULL;

static inline long time_diff(const struct timespec *start,
	const struct timespec *stop);
static int init_barrier(void);
static struct threads *init_tls(void);
static struct accel_particle *init_particles(void);
static void *thread_main(void *args);
static int thread_init(unsigned id);
static void thread_deinit(unsigned id);
static int thread_step(struct thread_state *state, unsigned step, long *us);
static int build_step(unsigned step, float radius, long *us);
static void sync_tree_particles(struct accel_particle tree_particles[],
	const struct particle_slice *slice);

static const size_t kib		   = (size_t)1 << 10;
static const size_t mib		   = kib << 10;
static const size_t gib		   = mib << 10;
static const size_t arena_size = 4 * gib;

static inline bool
step_continue(unsigned step)
{
	return options.steps == 0 || step < options.steps;
}

#include <unistd.h>

int
main(int argc, char *argv[argc])
{
	int res;
	if ((res = options_parse(argc, argv)))
		return (res == BHE_EARLY_EXIT) ? 0 : res;

#ifdef RENDER
	if ((res = render_init()))
		return res;
#endif // RENDER

	// Initialize the global (shared) state.

	if (unlikely((res = init_barrier())))
		return ENOMEM;
	if (unlikely((res = arena_init(&arena, arena_size))))
		return res;
	if (unlikely((particles = init_particles()) == NULL))
		return ENOMEM;
	if (unlikely((tls = init_tls()) == NULL))
		return ENOMEM;

	// Spawn p - 1 additional worker threads.
	const unsigned pthreads = options.threads - 1;
	pthread_t *threads		= malloc(sizeof(pthread_t) * pthreads);
	if (unlikely(threads == NULL))
		return ENOMEM;

	unsigned t;
	for (t = 0; t < pthreads; t++) {
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		res = pthread_create(&threads[t], &attr, thread_main,
			(void *)(uintptr_t)(t + 1));
		pthread_attr_destroy(&attr);

		if (res)
			goto exit;
	}

	// Init the main thread state.
	if ((res = thread_init(0))) {
		atomic_store_explicit(&thread_errno, res, memory_order_release);
		goto exit;
	}

	if (options.verbose)
		fprintf(stderr, "begin simulation ...\n");
	else
		// Print only the CSV file header.
		fprintf(stdout, "step,build,simulate\n");

	struct thread_state *state = &tls->states[0];
	for (unsigned step = 0; step_continue(step); step++) {
		long build_us, step_us;
		if ((res = build_step(step, state->radius, &build_us)))
			goto exit;
		if ((res = thread_step(state, step, &step_us)))
			goto exit;

		if (options.verbose)
			fprintf(stderr,
				"step t = %u:\n"
				"\tbuilt tree in: %ld us\n"
				"\tsimulation in: %ld us\n",
				step, build_us, step_us);
		else
			; // fprintf(stdout, "%u,%ld,%ld\n", step, build_us, step_us);

		// Recalculate the radius for the next iteration step.
		//
		// All other threads will either wait in barrier or exit, so it is safe
		// to iterate and update each thread's radius.
		float max_radius = 0.0;
		for (unsigned t = 0; t < options.threads; t++) {
			const float radius = tls->states[t].radius;
			if (radius > max_radius)
				max_radius = radius;
		}

		for (unsigned t = 0; t < options.threads; t++)
			tls->states[t].radius = max_radius;

#ifdef RENDER
		if (render_scene(particles, tls->states[0].radius))
			goto exit;
		usleep(500000);
#endif // RENDER
	}

exit:
	for (unsigned i = 0; i < t; i++) {
		void *thread_res;
		pthread_join(threads[i], &thread_res);

		int tres = (int)((uintptr_t)thread_res);
		if (tres > 0)
			fprintf(stderr, "Error in joined thread %d: %s\n", i,
				strerror(tres));
		thread_deinit(i);
	}

	free(threads);
	free(tls);
	free(particles);
	arena_deinit(&arena);

#ifdef RENDER
	render_deinit();
#endif // RENDER

	return (res != BHE_EARLY_EXIT) ? res : 0;
}

int
arena_init(struct arena *arena, size_t size)
{
	arena->memory = mmap(NULL, size, PROT_READ | PROT_WRITE,
		MAP_PRIVATE | MAP_ANON, -1, 0);
	if (unlikely(arena->memory == MAP_FAILED))
		return errno;

	arena->curr = arena->memory;
	arena->end	= (void *)((uintptr_t)arena->memory + size);

	return 0;
}

void
arena_deinit(struct arena *arena)
{
	const size_t size = (uintptr_t)arena->end - (uintptr_t)arena->curr;
	if (munmap(arena->memory, size))
		fprintf(stderr, "Failed to unmap arena: %s\n", strerror(errno));
}

static inline long
time_diff(const struct timespec *start, const struct timespec *stop)
{
	return (stop->tv_sec - start->tv_sec) * 1e6
		+ ((stop->tv_nsec - start->tv_nsec) / 1e3);
}

static int
init_barrier(void)
{
	pthread_barrierattr_t attr;
	int res;

	if ((res = pthread_barrierattr_init(&attr)))
		return res;
	if ((res = pthread_barrier_init(&barrier, &attr, options.threads)))
		return res;

	return 0;
}

static struct threads *
init_tls(void)
{
	struct threads *tls = malloc(sizeof(struct threads)
		+ (sizeof(struct thread_state) * options.threads));
	if (unlikely(tls == NULL))
		return NULL;

	tls->len = options.threads;

	return tls;
}

static struct accel_particle *
init_particles(void)
{
	if (options.seed != 0)
		srandom(options.seed);

	struct accel_particle *particles
		= malloc(sizeof(struct accel_particle) * options.particles);
	if (unlikely(particles == NULL))
		return NULL;

	verbose_printf("randomizing %zu particles within radius %.3f.\n",
		options.particles, options.radius);
	randomize_particles(particles, options.radius);
	verbose_printf("particle randomization complete.\n");

	return particles;
}

static void *
thread_main(void *args)
{
	int res;

	const unsigned id = (unsigned)((uintptr_t)args);
	if ((res = thread_init(id))) {
		atomic_store_explicit(&thread_errno, res, memory_order_release);
		return (void *)((uintptr_t)res);
	}

	struct thread_state *state = &tls->states[id];
	for (unsigned step = 0; step_continue(step); step++)
		if ((res = thread_step(state, step, NULL)))
			return (void *)((uintptr_t)res);

	return NULL;
}

static int
thread_init(unsigned id)
{
	struct thread_state *state = &tls->states[id];

	if (id == 0)
		state->particles = NULL;
	else {
		const size_t size = sizeof(struct accel_particle) * options.particles;
		if (unlikely((state->particles = malloc(size)) == NULL))
			return ENOMEM;
	}

	const size_t len   = options.particles / options.threads;
	const size_t rem   = options.particles % options.threads;
	const size_t start = (size_t)id * len;

	state->id	 = id;
	state->slice = (struct particle_slice) {
		.offset = start,
		.len	= (id == options.threads - 1) ? len + rem : len,
		.from	= (id == 0) ? particles : &state->particles[start],
	};
	state->radius = options.radius;

	if (id != 0)
		sync_tree_particles(state->particles, NULL);

	return 0;
}

static void
thread_deinit(unsigned id)
{
	struct thread_state *state = &tls->states[id];
	free(state->particles);
}

static int
thread_step(struct thread_state *state, unsigned step, long *us)
{
	struct timespec start, stop;

	pthread_barrier_wait(&barrier);

	if (atomic_load_explicit(&thread_errno, memory_order_acquire))
		return BHE_EARLY_EXIT;

	if (state->id == 0)
		clock_gettime(CLOCK_MONOTONIC, &start);

	float radius = particle_tree_simulate(&tree, &state->slice);
	memcpy(&particles[state->slice.offset], state->slice.from,
		sizeof(struct accel_particle) * state->slice.len);
	state->radius = radius;

	if (state->id == 0) {
		clock_gettime(CLOCK_MONOTONIC, &stop);
		*us = time_diff(&start, &stop);
	}

	// Wait for all threads to complete the current simulation step and
	// propagate their results, before synchronizing the global particle slice
	// with the thread's local one.
	pthread_barrier_wait(&barrier);

	if (state->id != 0)
		sync_tree_particles(state->particles, &state->slice);

	return 0;
}

static int
build_step(unsigned step, float radius, long *us)
{
	struct timespec start, stop;
	int res;

	clock_gettime(CLOCK_MONOTONIC, &start);
	if (options.optimize && step % 10 == 0)
		sort_particles(particles);

	res = particle_tree_build(&tree, particles, radius, &arena);
	if (unlikely(res)) {
		atomic_store_explicit(&thread_errno, res, memory_order_release);
		pthread_barrier_wait(&barrier);
		return res;
	}

	clock_gettime(CLOCK_MONOTONIC, &stop);
	*us = time_diff(&start, &stop);

	return 0;
}

static void
sync_tree_particles(struct accel_particle tree_particles[],
	const struct particle_slice *slice)
{
	if (slice == NULL) {
		memcpy(tree_particles, particles,
			sizeof(struct accel_particle) * options.particles);
		return;
	}

	// Copy everything before the slice over into the local particle slice.
	memcpy(&tree_particles[0], &particles[0],
		sizeof(struct accel_particle) * slice->offset);
	// Copy everything after the slice over into the local particle slice.
	const size_t after = slice->offset + slice->len;
	const size_t len   = options.particles - after;
	memcpy(&tree_particles[after], &particles[after],
		sizeof(struct accel_particle) * len);
}
