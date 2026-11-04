#ifndef BARNES_HUT_THREAD_H
#define BARNES_HUT_THREAD_H

#include <stddef.h>
#include <threads.h>

#include "barnes-hut/def.h"
#include "barnes-hut/types.h"
#include "barnes-hut/physics.h"

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
	atomic int	state;
	cnd_t		cnd;
	mtx_t		mtx;
};

// fixme: need one wake/wait counter for the pool
// fixme: need N wake/wait state for each worker

struct thread_pool_worker {
	struct thread_pool			*pool;
	thrd_t						thread;
	bool						is_ready;

	struct thread_notify		state;

	struct particle				*particles;
	size_t						particle_count;
	size_t						task_size;
	size_t						sort_task_task;
	error_t						error;
	scalar_t					len;
};

struct thread_pool {
	struct thread_pool_worker	main_worker;

	struct {
		struct thread_notify	running_count;
	} aligned(64);

	struct {
		volatile size_t			work_idx;
	} aligned(64);

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
