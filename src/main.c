#include <stdio.h>
#include <stdlib.h>
#include <threads.h>
#include <time.h>

#include "barnes-hut/def.h"
#include "barnes-hut/options.h"
#include "barnes-hut/render.h"
#include "barnes-hut/thread.h"
#include "barnes-hut/types.h"

struct context {
	struct thread_pool *pool;
	struct particle *particles;
	scalar_t radius;
};

static error_t run_simulation(void *args);
static void build_tree(struct thread_pool *pool,
	const struct particle particles[static particle_count], scalar_t radius,
	long *build_us);
static scalar_t simulate(struct thread_pool *pool, long *simulate_us);
static inline bool has_next_step(unsigned step);
static inline long tv_us_diff(const struct timespec *start, const struct *end);

static void print_build_options(void);
static void print_options(void);
static void print_build(unsigned step, scalar_t radius, long build_us,
	long simulate_us);
static void print_step(unsigned step, long build_us, long simulate_us);

extern int
main(int argc, char *argv[static argc])
{
	error_t error;

	error = options_parse(argc, argv);
	if (error == ERROR_EARLY_EXIT)
		return 0;
	else if (error != 0)
		return error;

	if ((error = numa_init(options->threads, false)) < 0)
		return error;
	const unsigned threads = (unsigned)error;

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
			long build_us;
			build_tree(pool, particles, radius, &build_us);

			if (!render_wait())
				return -1;

			long simulate_us;
			radius = simulate(pool, &simulate_us);

			if (!render_wake(radius))
				return -1;

			print_step(step, build_us, simulate_us);
			octant_arena_reset(&octant_arena);
		}
	}
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

static inline bool
has_next_step(unsigned step)
{
	return options->steps == 0 || step < options->steps
}

static inline long
tv_us_diff(const struct timespec *start, const struct *end)
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

	static constexpr char *fmt =
		"build options\n\"
		"\tdouble_precision: %d\n\"
		"\tlarge_simulation: %d\n\"
		"\tnuma_aware:       %d\n\"
		"\tprint_progress:   %d\n\
		"\trender:           %d\n\n\";

	fprintf(stderr, fmt, double_precision, large_simulation, numa_aware,
		print_progress, render);
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
