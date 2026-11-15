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

extern int
main(int argc, char *argv[static argc])
{
	error_t error;

	if ((error = options_parse(argc, argv)))
		return error;
	// fixme: early exit

	if ((error = numa_init(options->threads, false)) < 0)
		return error;
	const unsigned threads = (unsigned)error;

	struct thread_pool *pool = thread_pool_init(threads);
	struct context ctx = {
		.pool = pool,
		.particles = particles,
		.radius = options->radius,
	};

#ifdef OPTION_RENDER
	if ((error = render_main())) {
		// join simulation thread
		return error;
	}

#else

#endif /* OPTION_RENDER */
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
		struct timespec start, end;

		clock_gettime(CLOCK_MONOTONIC, &start);
		thread_pool_build(pool, particles, radius);
		clock_gettime(CLOCK_MONOTONIC, &end);
		const long build_us = tv_us_diff(&start, &end);

		if (!render_wait())
			return -1;

		clock_gettime(CLOCK_MONOTONIC, &start);
		radius = thread_pool_simulate(pool);
		clock_gettime(CLOCK_MONOTONIC, &end);
		const long simulate_us = tv_us_diff(&start, &end);

		if (!render_wake(radius))
			return -1;

		print_step(step, build_us, simulate_us);
		octant_arena_reset(&octant_arena);
	}
}

static void print_options(void);
static void print_build_options(void);
static void print_build(unsigned step, scalar_t radius, long build_us, long simulate_us);
static void print_step(unsigned step, long build_us, long simulate_us);

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
