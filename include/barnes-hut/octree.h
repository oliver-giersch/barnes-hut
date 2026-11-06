#ifndef BARNES_HUT_OCTREE_H
#define BARNES_HUT_OCTREE_H

#include <stdatomic.h>
#include <stddef.h>
#include <stdint.h>

#include "barnes-hut/def.h"
#include "barnes-hut/numa.h"
#include "barnes-hut/physics.h"
#include "barnes-hut/types.h"

typedef uint32_t octant_id_t;

/* The size of the octant arena. */
extern size_t octant_count;

static constexpr size_t octant_children = 8;
static constexpr octant_id_t invalid_id = UINT32_MAX;

enum octant_flags_t : uint32_t {
	OCTANT_IS_LEAF = 0,
	OCTANT_IS_ABSORBING = 1,
	OCTANT_IS_NODE = 2,
};
struct octant {
	struct point_mass center;
	struct vec3 corner;
	_Atomic enum octant_flags_t is_node;
	octant_id_t children[octant_children];
};

extern void octant_init(struct octant *octant, const struct point_mass *center,
	scalar_t x, scalar_t y, scalar_t z);
extern void octant_insert_threadsafe(struct octant *octant,
	struct point_mass *pm, scalar_t len, struct octant_cache *cache);

struct octant_arena_node {
	octant_id_t start, end;
	octant_id_t idx;
} numa_aligned;

struct octant_arena {
	struct octant *memory;
	unsigned node_count;
	struct octant_arena_node nodes[numa_max_nodes];
};

/* The global octant arena. */
extern struct octant_arena octant_arena;

/* ... */
static constexpr size_t octant_cache_len = 32;

struct octant_cache {
	size_t idx;
	octant_id_t array[octant_cache_len];
};

static inline void octant_arena_reset(struct octant_arena *arena);
static inline octant_id_t octant_arena_alloc(struct octant_arena *arena);
static inline octant_id_t octant_arena_alloc_n_threadsafe(
	struct octant_arena *arena, unsigned n);
static inline struct octant *octant_arena_resolve(struct octant_arena *arena,
	octant_id_t id);

static inline void
octant_arena_reset(struct octant_arena *arena)
{
	for (unsigned i = 0; i < arena->node_count; i++) {
		arena->nodes[i].idx = arena->nodes[i].start;
	}
}

static inline octant_id_t
octant_arena_alloc(struct octant_arena *arena)
{
	struct octant_arena_node *node = &arena->nodes[numa_node_id];
	if (node->idx >= node->end)
		return invalid_id;

	octant_id_t id = node->idx;
	node->idx += 1;

	return id;
}

static inline octant_id_t
octant_arena_alloc_n_threadsafe(struct octant_arena *arena, unsigned n)
{
	struct octant_arena_node *node = &arena->nodes[numa_node_id];
	octant_id_t atomic *ptr = (octant_id_t atomic *)&node->idx;

	if (atomic_load_explicit(ptr, memory_order_relaxed) >= node->end)
		return invalid_id;

	octant_id_t id = atomic_fetch_add_explicit(ptr, n, memory_order_relaxed);
	octant_id_t end = id + n;

	if (end >= node->end)
		return invalid_id;

	return id;
}

static inline struct octant *
octant_arena_resolve(struct octant_arena *arena, octant_id_t id)
{
	// assert(id != invalid_octant && id < ocant_count, "...");
	return &arena->memory[id];
}

#endif /* BARNES_HUT_OCTREE_H */
