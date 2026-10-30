#ifndef BARNES_HUT_THREAD_H
#define BARNES_HUT_THREAD_H

#include <stddef.h>
#include <threads.h>

#include "barnes-hut/types.h"

enum thread_state : int {
	THREAD_IDLE			= 0,
	THREAD_EXIT			= 1,
	THREAD_RANDOMIZE	= 2,
	THREAD_SORT			= 3,
	THREAD_BUILD		= 4,
	THREAD_SIMULATE		= 5,
};

struct thread_pool;

struct thread_notify {
	cnd_t	cnd;
	mtx_t	mtx;
};

// fixme: need one wake/wait counter for the pool
// fixme: need N wake/wait state for each worker

struct thread_pool_worker {
	struct thread_pool			*pool;
	bool						is_ready;

	struct thread_notify		notify;
	volatile enum thread_state	state;

	struct particle				*particles;
	size_t						particle_count;
	size_t						task_size;
	size_t						sort_task_task;
	error_t						error;
};

struct thread_pool {
	struct thread_pool_worker	main_worker;

	struct thread_notify		notify;
	size_t						running_count;

	volatile size_t				work_idx;

	size_t						worker_count;
	struct thread_pool_worker	*workers[];
};

extern struct thread_pool *thread_pool_init(size_t threads);
extern void thread_pool_deinit(struct thread_pool *pool);
extern void thread_pool_randomize(struct thread_pool *pool,
	struct particle particles[], scalar_t radius);
extern scalar_t thread_pool_build(struct thread_pool *pool,
	struct particle particles[], scalar_t radius);
extern scalar_t thread_pool_simulate(struct thread_pool *pool);

#endif /* BARNES_HUT_THREAD_H */

