#ifndef BARNES_HUT_OCTREE_H
#define BARNES_HUT_OCTREE_H

#include <stddef.h>
#include <stdint.h>

#include "barnes-hut/physics.h"

typedef uint32_t octant_id_t;

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
	enum octant_flags_t is_node;
	octant_id_t children[octant_children];
};

struct octant_arena {
	struct octant *memory;
	octant_id_t idx;
	octant_id_t len;
};

/* The global octant arena. */
extern struct octant_arena arena;

/* ... */
static constexpr size_t octant_cache_len = 32;

struct octant_cache {
	octant_id_t array[octant_cache_len];
	size_t idx;
};

static inline void octant_arena_reset(struct octant_arena *arena);
static inline octant_id_t octant_arena_alloc(struct octant_arena *arena);
static inline struct octant *octant_arena_resolve(struct octant_arena *arena,
	octant_id_t id);

#endif /* BARNES_HUT_OCTREE_H */
