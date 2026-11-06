#ifndef BARNES_HUT_NUMA_H
#define BARNES_HUT_NUMA_H

#include <stddef.h>

#include "barnes-hut/types.h"

#ifdef OPTION_NUMA_AWARE
#include "barnes-hut/def.h"

#define numa_aligned aligned(cache_line_size)

static constexpr size_t numa_max_nodes = 8;

extern thread_local unsigned numa_node_id;

extern int numa_init(void);
extern int numa_thread_init(unsigned id);
extern void numa_destroy(void);
#else
#define numa_aligned

static constexpr size_t numa_max_nodes = 1;
static constexpr unsigned numa_node_id = 0;

static inline int numa_init(void)
{
	return (ERROR_SUCCESS);
}

static inline int
numa_thread_init(unsigned id)
{
	(void)id;
	return (ERROR_SUCCESS);
}

static inline void
numa_destroy(void) {}
#endif /* OPTION_NUMA_AWARE */
#endif /* BARNES_HUT_NUMA_H */
