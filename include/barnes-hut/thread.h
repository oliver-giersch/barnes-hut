#ifndef BARNES_HUT_THREAD_H
#define BARNES_HUT_THREAD_H

#include <stddef.h>
#include <threads.h>

#include "barnes-hut/def.h"
#include "barnes-hut/futex.h"
#include "barnes-hut/numa.h"
#include "barnes-hut/physics.h"
#include "barnes-hut/types.h"

enum thread_state : int {
	THREAD_IDLE = 0,
	THREAD_EXIT = 1,
	THREAD_RANDOMIZE = 2,
	THREAD_SORT = 3,
	THREAD_BUILD = 4,
	THREAD_RECENTER = 5,
	THREAD_SIMULATE = 6,
};

struct thread_pool;

struct thread_notify {
	atomic int state;
	cnd_t cnd;
	mtx_t mtx;
};

// fixme: need one wake/wait counter for the pool
// fixme: need N wake/wait state for each worker

struct thread_pool_worker {
	struct {
		struct thread_pool *pool;
		thrd_t thread;
		unsigned idx;
		bool is_ready;
	} attr_cache_line_aligned;

	struct {
		atomic_uint32_t state;
	} attr_cache_line_aligned;

	struct {
		struct particle *particles;
		struct octant *root;
		size_t particle_count;
		error_t error;
		scalar_t len;
	} attr_cache_line_aligned;
};

struct thread_pool {
	/* Thread pool constants. */
	struct {
		/* The global particle array. */
		struct particle *particles;
		/* ... */
		size_t step_task_size;
		/* ... */
		size_t worker_count;
	} attr_cache_line_aligned;

	struct {
		struct thread_pool_worker main_worker;
		atomic_uint32_t running_count;
	} attr_cache_line_aligned;

	struct {
		atomic size_t work_idx;
		size_t start, end;
	} attr_cache_line_aligned task_nodes[numa_max_nodes];

	struct thread_pool_worker *workers[];
};

extern struct thread_pool *thread_pool_init(
	struct particle particles[static particle_count], size_t threads);
extern void thread_pool_destroy(struct thread_pool *pool);
extern void thread_pool_randomize(struct thread_pool *pool,
	struct particle particles[static particle_count], scalar_t radius);
extern scalar_t thread_pool_build(struct thread_pool *pool,
	const struct particle particles[static particle_count], scalar_t radius);
extern scalar_t thread_pool_simulate(struct thread_pool *pool);

#endif /* BARNES_HUT_THREAD_H */
