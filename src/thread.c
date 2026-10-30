#include "barnes-hut/thread.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <threads.h>

#include "barnes-hut/physics.h"

static constexpr size_t thread_pool_alignment = 4096;
static constexpr size_t pool_mask = thread_pool_alignment - 1;

static error_t thread_pool_worker_main(void *arg);

/* Notify */

typedef bool (*thread_notify_check_fn)(void *);
typedef void (*thread_notify_update_fn)(void *);

static error_t thread_notify_init(struct thread_notify *notify);
static void thread_notify_wait(struct thread_notify *notify,
	thread_notify_check_fn fn, void *ctx);
static void thread_notify_wake(struct thread_notify *notify,
	thread_notify_update_fn fn, void *ctx);
static void thread_notify_wake_one(struct thread_notify *notify,
	thread_notify_update_fn fn, void *ctx);

typedef struct {
	cnd_t cnd;
	mtx_t mtx;
	int state;
} futex_t;

static inline error_t
futex_init(futex_t *futex, int init)
{
	error_t error;

	if ((error = cnd_init(&futex->cnd)))
		return (error);
	if ((error = mtx_init(&futex->mtx, mtx_plain)))
		return (error);

	futex->state = init;
	return (0);
}

static inline void
futex_wait(futex_t *futex, int expected)
{
	mtx_lock(&futex->mtx);

	while (futex->state == expected)
		cnd_wait(&futex->cnd, &futex->mtx);

	mtx_unlock(&futex->mtx);
	return;
}

typedef int(update_fn_t)(int);

static inline void
futex_wake_one(futex_t *futex, update_fn_t fn)
{
	mtx_lock(&futex->mtx);
	futex->state = (fn)(futex->state);
	cnd_signal(&futex->cnd);
	mtx_unlock(&futex->mtx);
}

static inline void
futex_wake(futex_t *futex, update_fn_t fn)
{
	mtx_lock(&futex->mtx);
	futex->state = (fn)(futex->state);
	cnd_broadcast(&futex->cnd);
	mtx_unlock(&futex->mtx);
}

extern struct thread_pool *
thread_pool_init(size_t threads)
{
	size_t thread_pool_size = sizeof(struct thread_pool)
		+ (threads * sizeof(struct thread_pool_worker *));

	thread_pool_size = (thread_pool_size <= thread_pool_alignment)
		? thread_pool_alignment
		: (thread_pool_size / thread_pool_alignment) * thread_pool_alignment;

	struct thread_pool *pool
		= aligned_alloc(thread_pool_alignment, thread_pool_size);
	if (pool == nullptr)
		return nullptr;

	pool->main_worker = (struct thread_pool_worker) {
		.pool = pool,
		.state = THREAD_IDLE,
		.is_ready = false,
	},
	pool->worker_count = threads;
	pool->workers[0] = &pool->main_worker;

	return pool;
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

	for (int i = 0; i < pool->worker_count; i++)
		pool->worker[i].len = len;

	pool->work_idx = 1; // FIXME: atomic store?
	thread_pool_run(pool, THREAD_BUILD);

	for (int i = 0; i < pool->worker_count; i++)
		if (pool->workers[i].error != 0)
			return FP_NAN; // FIXME: return int/error

	// FIXME: (parallel?) update center?
}

static error_t
thread_pool_worker_main(void *arg)
{
	const uintptr_t ptr = (uintptr_t)arg;

	struct thread_pool *pool = (struct thread_pool *)(ptr & ~pool_mask);
	const int idx = (int)(ptr & pool_mask);

	struct thread_pool_worker worker = {
		.pool = pool,
		.state = THREAD_IDLE,
		.is_ready = false,
	};

	pool->workers[idx] = &worker;

	// atomic release-store
	// futex_wake_one

	// then loop forever, waiting for tasks
}

static error_t
thread_notify_init(struct thread_notify *notify)
{
	error_t error;

	if ((error = cnd_init(&notify->cnd)))
		return error;
	if ((error = mtx_init(&notify->mtx, mtx_plain)))
		return error;
	return 0;
}

static void
thread_notify_wait(struct thread_notify *notify, thread_notify_check_fn fn,
	void *ctx)
{
	mtx_lock(&notify->mtx);

	if ((fn)(ctx))
		cnd_wait(&notify->cnd, &notify->mtx);

	mtx_unlock(&notify->mtx);
}

static void thread_notify_wake(struct thread_notify *notify,
	thread_notify_update_fn fn, void *ctx);
static void thread_notify_wake_one(struct thread_notify *notify,
	thread_notify_update_fn fn, void *ctx);
