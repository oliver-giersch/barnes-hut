#include "barnes-hut/thread.h"

#include <math.h>
#include <stdatomic.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <threads.h>

#include "barnes-hut/def.h"
#include "barnes-hut/octree.h"
#include "barnes-hut/options.h"
#include "barnes-hut/physics.h"
#include "barnes-hut/types.h"

static constexpr size_t thread_pool_alignment = 4096;
static constexpr size_t pool_mask = thread_pool_alignment - 1;

static void thread_pool_run(struct thread_pool *pool, enum thread_state task);
static void thread_pool_wait_done(struct thread_pool *pool);
static void thread_pool_notify_done(struct thread_pool *pool);

struct batch {
	struct particle	*start;
	struct particle	*end;
};

static bool thread_pool_next_batch(struct thread_pool *pool,
	struct particle particles[], struct batch *batch);

/* Worker thread */

static error_t thread_pool_worker_main(void *arg);
static void thread_pool_worker_build(struct thread_pool_worker *worker,
	struct particle particles[], struct octant *root, scalar_t len);

static bool thread_pool_worker_next_batch(struct thread_pool_worker *worker,
	struct batch *batch);

/* Notify */

static error_t thread_notify_init(struct thread_notify *notify);
static void thread_notify_wait(struct thread_notify *notify, int expected);
static void thread_notify_wake_one(struct thread_notify *notify);
static void thread_notify_wake_all(struct thread_notify *notify);

/* Utility */

struct thread_ctx {
	struct thread_pool	*pool;
	unsigned			idx;
};

static inline void *thread_ctx_compose(struct thread_pool *pool, unsigned idx);
static inline struct thread_ctx thread_ctx_decompose(void *ctx);

/* Progress printing */

#if (OPTION_PRINT_PROGRESS == 1)
struct progress_bar {
	atomic size_t	*global_progress;
	size_t			local_progress;
	size_t			step;
};
#else
struct progress_bar {
	char _;
};
#endif /* OPTION_PRINT_PROGRESS */

static void progress_bar_print(struct progress_bar *bar,
	const struct thread_pool_worker *worker);

extern struct thread_pool *
thread_pool_init(size_t threads)
{
	error_t error;

	/* Calculate the (dynamic) thread pool memory size. */
	size_t thread_pool_size = sizeof(struct thread_pool)
		+ (threads * sizeof(struct thread_pool_worker *));

	thread_pool_size = (thread_pool_size <= thread_pool_alignment)
		? thread_pool_alignment
		: (thread_pool_size / thread_pool_alignment) * thread_pool_alignment;

	struct thread_pool *pool
		= aligned_alloc(thread_pool_alignment, thread_pool_size);
	if (pool == nullptr)
		return nullptr;

	/* Initialize the pool and main worker notify instances. */
	if ((error = thread_notify_init(&pool->running_count)))
		goto error_exit;
	if ((error = thread_notify_init(&pool->main_worker.state)))
		goto error_exit;

	/* Initialize the main thread worker and pool. */
	pool->worker_count = threads;
	pool->workers[0] = &pool->main_worker;
	pool->running_count.state = threads - 1;

	pool->main_worker.pool = pool;
	pool->main_worker.thread = thrd_current();
	pool->main_worker.is_ready = false;
	pool->main_worker.state.state = THREAD_IDLE;

	for (unsigned i = 1; i < threads; i++) {
		error = thrd_create(&pool->workers[i]->thread, thread_pool_worker_main,
			thread_ctx_compose(pool, i));
	}

	/* Wait until all worker threads have initialized their state. */
	thread_pool_wait_done(pool);

	return pool;

error_exit:
	free_aligned_sized(pool, thread_pool_alignment, thread_pool_size);
	return nullptr;
}

extern scalar_t
thread_pool_build(struct thread_pool *pool, struct particle particles[],
	scalar_t radius)
{
	const scalar_t corner = -radius;
	const scalar_t id = octant_arena_alloc(&octant_arena);
	const scalar_t root = octant_arena_resolve(&octant_arena, id);
	octant_init(root, &particles[0].pm, corner, corner, corner);

	const scalar_t len = radius * 2.0;

	for (unsigned i = 0; i < pool->worker_count; i++)
		pool->workers[i]->len = len;

	pool->work_idx = 1; // FIXME: atomic store?
	thread_pool_run(pool, THREAD_BUILD);

	for (unsigned i = 0; i < pool->worker_count; i++)
		if (pool->workers[i]->error != 0)
			return FP_NAN; // FIXME: return int/error

	// FIXME: (parallel?) update center?
}

static void
thread_pool_run(struct thread_pool *pool, enum thread_state task)
{
	// set run count
	// for each worker: set task and wake
}

static void
thread_pool_wait_done(struct thread_pool *pool)
{
	int running_count;

	while (true) {
		running_count = atomic_load_explicit(&pool->running_count.state,
			memory_order_relaxed);
		if (running_count == 0)
			return;

		thread_notify_wait(&pool->running_count, running_count);
	}
}

static void
thread_pool_notify_done(struct thread_pool *pool)
{
	const int running_count = atomic_fetch_sub_explicit(
		&pool->running_count.state, 1, memory_order_release) - 1;
	assert(running_count >= 0, "...");
	if (running_count == 0)
		thread_notify_wake_one(&pool->running_count);
}

static bool
thread_pool_next_batch(struct thread_pool *pool, struct particle particles[],
	struct batch *batch)
{
	const size_t task_size = pool->step_task_size;
	const size_t start = atomic_fetch_add_explicit(&pool->work_idx,
		task_size, memory_order_relaxed);

	if (start > particle_count)
		return false;

	const size_t end = (start + task_size < particle_count)
		? start + task_size
		: particle_count;

	*batch = (struct batch){
		.start = particles + start,
		.end = particles + end,
	};

	return true;
}

/* Worker thread */

static error_t
thread_pool_worker_main(void *arg)
{
	struct thread_ctx ctx = thread_ctx_decompose(arg);
	struct thread_pool *pool = ctx.pool;
	unsigned idx = ctx.idx;

	struct thread_pool_worker worker;

	error_t error = thread_notify_init(&worker.state);
	if (error)
		goto out;

	worker.pool = pool;
	worker.state.state = THREAD_IDLE;
	worker.is_ready = false;

	/* Register the initialized worker with the pool. */
	pool->workers[idx] = &worker;
	thread_pool_notify_done(pool);


	// atomic release-store
	// futex_wake_one

	// then loop forever, waiting for tasks
	while (true) {
		const int state = atomic_load_explicit(&worker.state.state,
			memory_order_acquire);
		if (state == THREAD_IDLE) {
			thread_notify_wait(&worker.state, THREAD_IDLE);
			continue;
		} else if (state == THREAD_EXIT)
			break;

		assert(worker.is_ready, "");
		if (state == THREAD_RANDOMIZE)
			;
		else if (state == THREAD_SORT)
			;
		else if (state == THREAD_BUILD)
			// publish the "local" extent!
			thread_pool_worker_build(&worker, worker.particles, worker.root,
				worker.len);
		else if (state == THREAD_SIMULATE)
			;
		else
			__builtin_unreachable();


	}

out:
	thread_pool_notify_done(pool);
	return (error);
}

static void
thread_pool_worker_build(struct thread_pool_worker *worker,
	struct particle particles[], struct octant *root, scalar_t len)
{
	struct batch batch;
	struct particle *particle;

	while (thread_pool_next_batch(worker->pool, particles, &batch)) {
		for (particle = batch.start; particle < batch.end; particle++) {
			octant_insert_non_blocking(root, &particle->pm, len, nullptr);
			progress_print(&progress);
		}
	}
}

static bool
thread_pool_worker_next_batch(struct thread_pool_worker *worker,
	struct batch *batch)
{
	assert(worker->is_ready, "");

	const size_t task_size = worker->task_size;
	const size_t start = atomic_fetch_add_explicit(&worker->pool->work_idx,
		task_size, memory_order_relaxed);

	if (start > particle_count)
		return false;

	const size_t end = (start + task_size < particle_count)
		? start + task_size
		: particle_count;

	*batch = (struct batch){
		.start = worker->particles + start,
		.end = worker->particles + end,
	};

	return true;
}

/* Notify */

static error_t
thread_notify_init(struct thread_notify *notify)
{
	error_t error;

	if ((error = cnd_init(&notify->cnd)))
		return error;

	if ((error = mtx_init(&notify->mtx, mtx_plain))) {
		cnd_destroy(&notify->cnd);
		return error;
	}

	return 0;
}

static void
thread_notify_wait(struct thread_notify *notify, int expected)
{
	if (atomic_load_explicit(&notify->state, memory_order_relaxed) != expected)
		return;

	mtx_lock(&notify->mtx);
	if (atomic_load_explicit(&notify->state, memory_order_relaxed) == expected)
		cnd_wait(&notify->cnd, &notify->mtx);
	mtx_unlock(&notify->mtx);
}

static void
thread_notify_wake_one(struct thread_notify *notify)
{
	mtx_lock(&notify->mtx);
	cnd_signal(&notify->cnd);
	mtx_unlock(&notify->mtx);
}

static void
thread_notify_wake_all(struct thread_notify *notify)
{
	mtx_lock(&notify->mtx);
	cnd_broadcast(&notify->cnd);
	mtx_unlock(&notify->mtx);
}

/* Utility */

static inline void
*thread_ctx_compose(struct thread_pool *pool, unsigned idx)
{
	//assert(idx < thread_pool_alignment);
	uintptr_t composed = (uintptr_t)pool;
	composed |= (uintptr_t)idx;

	return (void *)composed;
}

static inline struct thread_ctx
thread_ctx_decompose(void *ctx)
{
	uintptr_t addr = (uintptr_t)ctx;
	struct thread_ctx thread_ctx =  {
		.pool = (struct thread_pool *)(addr & ~pool_mask),
		.idx = (unsigned)(addr & pool_mask),
	};

	return thread_ctx;
}

/* Progress */

static void
progress_bar_print(struct progress_bar *bar,
	const struct thread_pool_worker *worker)
{
	if (!options->verbose)
		return;

#if (OPTION_PRINT_PROGRESS == 1)
	bar->local_progress += 1;
	if (bar->local_progress < bar->step)
		return;

	atomic_fetch_add_explicit(bar->global_progress, bar->local_progress,
		memory_order_relaxed);
	bar->local_progress = 0;

	if (worker->idx != 0)
		return;

	static constexpr char fill = '#'; // '█'
#else
	(void)bar;
	(void)worker;
#endif /* OPTION_PRINT_PROGRESS */
}
