#include "barnes-hut/thread.h"

#include <math.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <threads.h>

#include "barnes-hut/def.h"
#include "barnes-hut/octree.h"
#include "barnes-hut/physics.h"
#include "barnes-hut/types.h"

static constexpr size_t thread_pool_alignment = 4096;
static constexpr size_t pool_mask = thread_pool_alignment - 1;

static void thread_pool_run(struct thread_pool *pool, enum thread_state task);
static void thread_pool_wait_done(struct thread_pool *pool);
static void thread_pool_notify_done(struct thread_pool *pool);

static error_t thread_pool_worker_main(void *arg);

/* Notify */

static error_t thread_notify_init(struct thread_notify *notify);
static void thread_notify_wait(struct thread_notify *notify, int expected);
static void thread_notify_wake_one(struct thread_notify *notify);
static void thread_notify_wake_all(struct thread_notify *notify);

/* Progress printing */

struct progress {
	atomic size_t	*global_progress;
	size_t			local_progress;
};

/* Utility */

struct thread_ctx {
	struct thread_pool	*pool;
	unsigned			idx;
};

static inline void *thread_ctx_compose(struct thread_pool *pool, unsigned idx);
static inline struct thread_ctx thread_ctx_decompose(void *ctx);

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
	const auto corner = -radius;
	const auto id = octant_arena_alloc(&octant_arena);
	const auto root = octant_arena_resolve(&octant_arena, id);
	octant_init(root, &particles[0].pm, corner, corner, corner);

	const auto len = radius * 2.0;

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

static error_t
thread_pool_worker_main(void *arg)
{
	struct thread_ctx ctx = thread_ctx_decompose(arg);
	struct thread_pool *pool = ctx.pool;
	unsigned idx = ctx.idx;

	struct thread_pool_worker worker = {
		.pool = pool,
		.state = THREAD_IDLE,
		.is_ready = false,
	};

	worker.pool = pool;
	worker.is_ready = false;
	error_t error = thread_notify_init(&worker.state);


	pool->workers[idx] = &worker;
	thread_pool_notify_done(pool);


	// atomic release-store
	// futex_wake_one

	// then loop forever, waiting for tasks
	while (true) {
		const int state = atomic_load_explicit(&worker.state.state,
			memory_order_relaxed);
		if (state == THREAD_IDLE) {
			thread_notify_wait(&worker.state, THREAD_IDLE);
			continue;
		} else if (state == THREAD_EXIT)
			break;

	}
}

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
