#define _XOPEN_SOURCE 700

#include <stddef.h>
#include <stdlib.h>

#include <sys/mman.h>

#include <pthread.h>
#include <string.h>

#include "arena.h"
#include "common.h"
#include "options.h"
#include "phys.h"

struct thread_state {
	unsigned id;
	struct particle_tree *tree;
	struct particle_slice slice;
	struct arena arena;
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

// The global thread synchronization barrier.
static pthread_barrier_t barrier;
// The globally shared and synchronized region of all simulated particles.
static struct moving_particles *particles;
// The TLS holding the state of all threads.
static struct threads {
	unsigned len;
	struct thread_state states[];
} *tls = NULL;

static inline float
randomf(void)
{
	return (float)random() / (float)RAND_MAX;
}

static int init_barrier();
static struct threads *init_tls(void);
static struct moving_particle *init_particles(void);
static int thread_init(unsigned id);
static int thread_step(struct thread_state *state, unsigned step);
static void *thread_main(void *args);

static inline void sync_particles(const struct particle_slice *slice);
static inline void sync_slice(const struct particle_slice *slice);

int
main(int argc, char *argv[argc])
{
	int res;
	if ((res = options_parse(argc, argv)))
		return res;

	// init barrier: pthread_barrier_init(&barrier, &attr, options.threads);
	if ((res = init_barrier()))
		return -1;
	if ((particles = init_particles()) == NULL)
		return -1;
	if ((tls = init_tls()) == NULL)
		return -1;

	// spawn p - 1 additional worker threads
	const unsigned pthreads = options.threads - 1;
	pthread_t *threads		= malloc(sizeof(pthread_t) * pthreads);
	if (threads == NULL)
		return -1;

	unsigned t;
	for (t = 0; t < pthreads; t++) {
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		res = pthread_create(&threads[t], &attr, thread_main, (void *)(t + 1));
		pthread_attr_destroy(&attr);

		if (res)
			goto kill;
	}

	if ((res = thread_init(0)))
		goto kill;

	struct thread_state *state = &tls->states[0];
	for (unsigned step = 0; step < options.steps; step++) {
		if ((res = thread_step(state, step)))
			goto kill;

		// All other threads will either wait in barrier or exit, so it is safe
		// to iterate and update each thread's radius.
		float max_radius = 0.0;
		for (unsigned t = 0; t < options.threads; t++) {
			const float radius = tls->states[t].radius;
			if (radius > radius)
				max_radius = radius;
		}

		for (unsigned t = 0; t < options.threads; t++)
			tls->states[t].radius = max_radius;
	}

	return 0;
kill:
	return res;
}

static void *
thread_main(void *args)
{
	const unsigned id = (unsigned)((uintptr_t)args);
	if (thread_init(id))
		return NULL; // TODO: more meaningful value

	struct thread_state *state = &tls->states[id];
	for (unsigned step = 0; step < options.steps; step++)
		if (thread_step(state, step))
			return NULL; // TODO: more meaningful value

	return NULL;
}

static int
thread_step(struct thread_state *state, unsigned step)
{
	int res;

	pthread_barrier_wait(&barrier);
	if (options.optimize && step % 10 == 0)
		particle_tree_sort(&state->tree);
	const float radius = state->radius;
	if ((res = particle_tree_build(&state->tree, radius, &state->arena)))
		return res;

	sync_particles(&state->slice);
	pthread_barrier_wait(&barrier);
	sync_slice(&state->slice);

	return 0;
}

static inline void
sync_particles(const struct particle_slice *slice)
{
	memcpy(&particles[slice->offset], slice->from,
		sizeof(struct moving_particle) * slice->len);
}

static inline void
sync_thread_slice(struct thread_state *state)
{
	const struct particle_slice *slice = &state->slice;
	memcpy(state->tree->);
}