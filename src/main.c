#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <threads.h>
#include <time.h>

#include <sys/mman.h>

#include "barnes-hut/def.h"
#include "barnes-hut/numa.h"
#include "barnes-hut/octree.h"
#include "barnes-hut/options.h"
#include "barnes-hut/physics.h"
#include "barnes-hut/render.h"
#include "barnes-hut/thread.h"
#include "barnes-hut/types.h"

struct context {
	struct thread_pool *pool;
	struct particle *particles;
	scalar_t radius;
};

/* The virtual memory range base start (16 TiB). */
static constexpr uintptr_t vm_range_start = 0x0000100000000000ULL;
/* The virtual memory range end start (65536 TiB). */
static constexpr uintptr_t vm_range_end = 0x0000800000000000ULL;

static error_t run_simulation(void *args);
static void sort_particles(struct thread_pool *pool,
	struct particle particles[static particle_count], long *sort_us);
static void build_tree(struct thread_pool *pool,
	const struct particle particles[static particle_count], scalar_t radius,
	long *build_us);
static scalar_t simulate(struct thread_pool *pool, long *simulate_us);
/* Allocates `size` bytes in an exclusive 16TiB range of virtual memory. */
static void *mmap_memory(size_t size, error_t *error);

static void print_build_options(void);
static void print_options(void);
static inline void print_header(unsigned step);
static inline void print_sort(long sort_us);
static inline void print_build(long build_us);
static inline void print_step(unsigned step, scalar_t radius, long sort_us,
	long build_us, long simulate_us);
printf_like_n(3, 4) static char *bump_printf(char **ptr, const char *end,
	const char *fmt, ...);

static inline size_t universe_size(size_t particles, unsigned threads);
static inline bool has_next_step(unsigned step);
static inline long tv_us_diff(const struct timespec *start,
	const struct timespec *end);
static inline double us2s(long us);

extern int
main(int argc, char *argv[static argc])
{
	error_t error;

	/* Parse program options in argv. */
	if ((error = options_parse(argc, argv)) == ERROR_EXIT_EARLY)
		return 0;
	else if (error != 0)
		return error;

	print_build_options();
	print_options();

	/* Determine NUMA topology (if enabled). */
	if ((error = numa_init(options->threads, NUMA_ROUND_ROBIN, false)) < 0)
		return error;
	const unsigned threads = (unsigned)error;

	/* Allocate (mmap) particle memory. */
	particle_count = options->particles;
	struct particle *particles = mmap_memory(particle_count, &error);
	if (particles == nullptr)
		return error;

	/* Allocate (mmap) octant memory. */
	const size_t arena_size = universe_size(particle_count, threads);
	struct octant *octants = mmap_memory(arena_size, &error);
	if (octants == nullptr)
		return error;

	/* Initialize the global octant bump allocator. */
	octant_count = arena_size;
	octant_arena_init(&octant_arena, octants);

	/* Initialize the thread pool. */
	struct thread_pool *pool = thread_pool_init(threads);
	if (pool == nullptr)
		goto destroy_numa;

	struct context ctx = {
		.pool = pool,
		.particles = particles,
		.radius = options->radius,
	};

#ifdef OPTION_RENDER
	thrd_t thread;
	if ((error = thrd_create(&thread, run_simulation, (void *)&ctx)))
		goto destroy_pool;

	if ((error = render_main())) {
		thrd_join(thread, &error);
		goto destroy_pool;
	}
#else
	if ((error = run_simulation((void *)&ctx)))
		goto destroy_pool;
#endif /* OPTION_RENDER */

destroy_pool:
	thread_pool_destroy(pool);
destroy_numa:
	numa_destroy();

	return error;
}

static error_t
run_simulation(void *args)
{
	struct context *ctx;

	ctx = (struct context *)args;

	struct thread_pool *pool = ctx->pool;
	struct particle *particles = ctx->particles;
	scalar_t radius = ctx->radius;

	unsigned step = 0;
	while (has_next_step(step)) {
		for (unsigned step = 0; has_next_step(step); step++) {
			long sort_us = -1;

			print_header(step);
			if (options->sort_steps != 0 && step % options->sort_steps == 0)
				sort_particles(pool, particles, &sort_us);
			print_sort(sort_us);

			long build_us;
			build_tree(pool, particles, radius, &build_us);
			print_build(build_us);

			if (!render_wait())
				return -1;

			long simulate_us;
			radius = simulate(pool, &simulate_us);

			if (!render_wake(radius))
				return -1;

			print_step(step, build_us, simulate_us);
		}
	}
}

static void
sort_particles(struct thread_pool *pool,
	struct particle particles[static particle_count], long *sort_us)
{
	// TODO: implement, dumb qsort?
	// TODO: half-dumb qsort -> each thread pre-sorts slice, then sort again on 1 core.

	*sort_us = 0;
}

static void
build_tree(struct thread_pool *pool,
	const struct particle particles[static particle_count], scalar_t radius,
	long *build_us)
{
	struct timespec start, end;

	clock_gettime(CLOCK_MONOTONIC, &start);
	thread_pool_build(pool, particles, radius);
	clock_gettime(CLOCK_MONOTONIC, &end);

	*build_us = tv_us_diff(&start, &end);
}

static scalar_t
simulate(struct thread_pool *pool, long *simulate_us)
{
	struct timespec start, end;

	clock_gettime(CLOCK_MONOTONIC, &start);
	const scalar_t radius = thread_pool_simulate(pool);
	clock_gettime(CLOCK_MONOTONIC, &end);

	*simulate_us = tv_us_diff(&start, &end);
	return radius;
}

static void *
mmap_memory(size_t size, error_t *error)
{
	static uintptr_t next_addr = vm_range_start;
	uintptr_t addr;

	*error = ENOMEM;
	for (addr = next_addr; addr <= vm_range_end; addr += vm_range_start) {
		void *ptr = mmap((void *)addr, size, PROT_READ | PROT_WRITE,
			MAP_PRIVATE | MAP_ANON | MAP_FIXED, -1, 0);

		if (ptr != MAP_FAILED) {
			next_addr = addr + vm_range_start;
			*error = 0;
			return ptr;
		}

		*error = errno;
	}

	return nullptr;
}

static inline size_t
universe_size(size_t particles, unsigned threads)
{
	const size_t min_size = threads * octant_cache_len;
	const size_t max_size = octant_arena_max_size;
	const size_t size = (particles < particles * 3)
		? particles * 3
		: max_size;

	if (size < min_size)
		return min_size;
	if (size > max_size)
		return max_size;
	return size;
}

static inline bool
has_next_step(unsigned step)
{
	return options->steps == 0 || step < options->steps;
}

static inline long
tv_us_diff(const struct timespec *start, const struct timespec *end)
{
	const long delta_sec = (long)(end->tv_sec - start->tv_sec);
	const long delta_nsec = (long)(end->tv_nsec - start->tv_nsec);
	return (delta_sec * (long)1e6) + (delta_nsec * (long)1e3);
}

static void
print_build_options(void)
{
	bool double_precision = false;
	bool large_simulation = false;
	bool numa_aware = false;
	bool print_progress = false;
	bool render = false;

	if (!options->verbose)
		return;

#ifdef OPTION_DOUBLE_PRECISION
	double_precision = true;
#endif /* OPTION_DOUBLE_PRECISION */

#ifdef OPTION_LARGE_SIMULATION
	large_simulation = true;
#endif /* OPTION_LARGE_SIMULATION */

#ifdef OPTION_NUMA_AWARE
	numa_aware = true;
#endif /* OPTION_NUMA_AWARE */

#ifdef OPTION_PRINT_PROGRESS
	print_progress = true;
#endif /* OPTION_PRINT_PROGRESS */

#ifdef OPTION_RENDER
	render = true;
#endif /* OPTION_RENDER */

	static constexpr char fmt[] =
		"build options\n"
		"\tdouble_precision: %d\n"
		"\tlarge_simulation: %d\n"
		"\tnuma_aware:       %d\n"
		"\tprint_progress:   %d\n"
		"\trender:           %d\n\n";

	fprintf(stderr, fmt, double_precision, large_simulation, numa_aware,
		print_progress, render);
}

static void
print_options(void)
{
	char buf[128];

	if (options->verbose)
		return;

	static constexpr char fmt[] =
		"options\n"
		"\tsteps:      %s\n"
		"\tparticles:  %zu\n"
		"\tmass:       %f\n"
		"\tradius:     %f\n"
		"\tthreads:    %s\n"
		"\tseed:       %d\n"
		"\tsort_steps: %s\n"
		"\tflat:       %s\n"
		"\tverbose:    %s\n"
		"\ttheta:      %f\n"
		"\tdt:         %f\n"
		"\tfile:       %s\n"
		"\tdump:       %s\n\n";

	static const char *booleans[] = {
		[false] = "false",
		[true] = "true",
	};

	char *ptr = buf;
	const char *steps = (options->steps > 0)
		? bump_printf(&ptr, buf + sizeof(buf), "%u", options->steps)
		: bump_printf(&ptr, buf + sizeof(buf), "%s", "unlimited");
	const char *threads = (options->threads > 0)
		? bump_printf(&ptr, buf + sizeof(buf), "%u", options->threads)
		: bump_printf(&ptr, buf + sizeof(buf), "%s", "maximum");
	const char *sort_steps = (options->sort_steps > 0)
		? bump_printf(&ptr, buf + sizeof(buf), "%u", options->sort_steps)
		: bump_printf(&ptr, buf + sizeof(buf), "%s", "disabled");

	fprintf(stderr, fmt, steps, options->particles, options->mass,
		options->radius, threads, options->seed, sort_steps,
		booleans[options->flat], booleans[options->verbose], options->theta,
		options->dt, options->file, booleans[options->dump]);
}

static void
print_header(unsigned step)
{
	if (options->verbose)
		return;

	fprintf(stderr, "step t = %u\n", step);
}

static inline void
print_sort(long sort_us)
{
	if (!options->verbose)
		return;

	if (sort_us < 0)
		fprintf(stderr, "\tsort:\t%15s %15s\n", "N/A", "N/A");
	else
		fprintf(stderr, "\tsort:\t%12ld us %8.2f s\n", sort_us, us2s(sort_us));
}

static inline void
print_build(long build_us)
{
	if (!options->verbose)
		return;

	fprintf(stderr, "\tbuild:\t%12ld us %8.2f s\n", build_us, us2s(build_us));
}

static inline void
print_step(unsigned step, scalar_t radius, long sort_us, long build_us,
	long simulate_us)
{
	fprintf(stdout, "%u,%ld,%ld,%ld\n", step, sort_us, build_us, simulat_us);
}

printf_like_n(3, 4)
static inline char *bump_printf(char **ptr, const char *end,
	const char *fmt, ...)
{
	if (*ptr == nullptr)
		return nullptr;

	const size_t room = (size_t)(end - *ptr);
	if (room == 0) {
		*ptr = nullptr;
		return nullptr;
	}

	va_list ap;
	va_start(ap, fmt);
	int n = vsnprintf(*ptr, room, fmt, ap);
	va_end(ap);

	if (n < 0 || (size_t)n >= room) {
		*ptr = nullptr;
		return nullptr;
	}

	char *str = *ptr;
	*ptr += n + 1;
	return str;
}

static inline double
us2s(long us)
{
	return (double)us / 1'000'000;
}

extern noreturn
__panic(const char *file, const char *function, int line, const char *condition,
	const char *message)
{
	if (message != NULL)
		fprintf(stderr, "%s:%d(%s): assert failed: %s (%s)", function, line,
			file, condition, message);
	else
		fprintf(stderr, "%s:%d(%s): assert failed: %s", function, line, file,
			condition);

	abort();
}
