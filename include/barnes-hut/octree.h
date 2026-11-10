#ifndef BARNES_HUT_OCTREE_H
#define BARNES_HUT_OCTREE_H

#include <stdatomic.h>
#include <stddef.h>
#include <stdint.h>

#include "barnes-hut/numa.h"
#include "barnes-hut/physics.h"
#include "barnes-hut/types.h"

#ifdef OPTION_HUGE_SIMULATION
typedef struct octant * octant_id_t;
static constexpr octant_id_t invalid_id = nullptr;
#else
typedef uint32_t octant_id_t;
static constexpr octant_id_t invalid_id = UINT32_MAX;
#endif /* OPTION_HUGE_SIMULATION */

/* The size of the particles array. */
extern size_t particle_count;
/* The size of the octant arena. */
extern size_t octant_count;

enum octant_flags_t : uint32_t {
	OCTANT_IS_LEAF = 0,
	OCTANT_IS_BUSY = 1,
	OCTANT_IS_NODE = 2,
};

static constexpr unsigned octant_children = 8;

struct octant {
	struct point_mass center;
	struct vec3 corner;
	_Atomic enum octant_flags_t flags;
	octant_id_t children[octant_children];
};

extern struct octant *octant_tree_init(
	const struct particle particles[static particle_count], scalar_t radius);
extern void octant_insert(struct octant *octant, const struct point_mass *pm,
	scalar_t len);
extern void octant_insert_threadsafe(struct octant *octant,
	const struct point_mass *pm, scalar_t len, struct octant_cache *cache);
extern struct vec3 octant_update_center(struct octant *octant);
extern void octant_update_accel(const struct octant *octant,
	const struct point_mass *pm, scalar_t len, scalar_t theta_threshold,
	struct vec3 *accel);

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
static constexpr unsigned octant_cache_len = 32;

struct octant_cache {
	size_t idx;
	octant_id_t array[octant_cache_len];
};
#endif /* BARNES_HUT_OCTREE_H */
