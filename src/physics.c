#include "barnes-hut/physics.h"

#include <math.h>
#include <stdatomic.h>
#include <string.h>
#include <tgmath.h>

#include "barnes-hut/def.h"
#include "barnes-hut/octree.h"
#include "barnes-hut/types.h"

static constexpr scalar_t EPSILON = 1e-6;
static const struct vec3 zero_vec = { 0.0, 0.0, 0.0 };

/* Octant */

struct octant_insertion;

static void octant_init(struct octant *octant, const struct point_mass *center,
	struct vec3 corner);
static void octant_insert_child(struct octant *octant,
	const struct point_mass *pm, scalar_t len);
static void octant_insert_leaf(struct octant *octant,
	const struct point_mass *pm, const struct octant_insertion *insertion);
static inline bool octant_is_leaf(const struct octant *octant);

static void octant_insert_child_threadsafe(struct octant *octant,
	const struct point_mass *pm, scalar_t child_len,
	struct octant_cache *cache);
static inline bool octant_is_leaf_threadsafe(const struct octant *octant);

static inline bool octant_absorbs(const struct octant *octant,
	const struct point_mass *pm, scalar_t child_len);


/* Octant insertion */

struct octant_insertion {
	struct vec3 corner;
	uint8_t children_idx;
};

static struct octant_insertion octant_insertion_prepare(
	const struct octant *octant, const struct point_mass *pm,
	scalar_t child_len);

/* Octant arena */

static inline void octant_arena_reset(struct octant_arena *arena);
static inline octant_id_t octant_arena_alloc(struct octant_arena *arena);
static inline octant_id_t octant_arena_alloc_n_threadsafe(
	struct octant_arena *arena, unsigned n);
static inline struct octant *octant_arena_resolve(struct octant_arena *arena,
	octant_id_t id);

/* Vector */

static inline struct vec3 vec3_scalar(scalar_t s);
static inline bool vec3_eql(const struct vec3 *restrict v,
	const struct vec3 *restrict u);
static inline scalar_t vec3_dist(const struct vec3 *restrict v,
	const struct vec3 *restrict u);
static inline scalar_t vec3_dist_sq(const struct vec3 *restrict v,
	const struct vec3 *restrict u);
static inline void vec3_assign_add(struct vec3 *restrict v,
	const struct vec3 *restrict u);
static inline void vec3_assign_sub(struct vec3 *restrict v,
	const struct vec3 *restrict u);
static inline void vec3_assign_mul(struct vec3 *v, scalar_t t);
static inline void vec3_assign_div(struct vec3 *v, scalar_t t);
static inline void vec3_assign_muladd(struct vec3 *restrict v,
	const struct vec3 *restrict u, scalar_t t);

static inline bool scalar_eql(scalar_t a, scalar_t b);
static inline scalar_t scalar_sq(scalar_t x);
static inline scalar_t scalar_random(void);
static inline scalar_t scalar_max(scalar_t a, scalar_t b);
static inline scalar_t scalar_sqrt(scalar_t x);
static inline scalar_t scalar_abs(scalar_t x);

extern struct octant *
octant_tree_init(const struct particle particles[static particle_count],
	scalar_t radius)
{
	const octant_id_t id = octant_arena_alloc(&octant_arena);
	if (unlikely(id == invalid_id))
		abort_oom();

	struct octant *root = octant_arena_resolve(&octant_arena, id);
	octant_init(root, &particles[0].pm, vec3_scalar(-radius));

	return (root);
}

extern void
octant_insert(struct octant *octant, const struct point_mass *pm, scalar_t len)
{
	const scalar_t child_len = len / 2.0;
	if (octant_is_leaf(octant)) {
		if (octant_absorbs(octant, pm, child_len)) {
			octant->center.mass += pm->mass;
			return;
		}

		const struct octant_insertion insertion = octant_insertion_prepare(
			octant, pm, child_len);
		octant_insert_leaf(octant, pm, &insertion);
		octant->flags = OCTANT_IS_NODE;
	}

	/* Recursively insert the particle. */
	octant_insert_child(octant, pm, child_len);
}

extern void
octant_insert_threadsafe(struct octant *octant, const struct point_mass *pm,
	scalar_t len, struct octant_cache *cache)
{
	const scalar_t child_len = len / 2.0;
}

extern struct vec3
octant_update_center(struct octant *octant)
{
	struct vec3 mass_moment;

	if (octant_is_leaf(octant)) {
		mass_moment = octant->center.pos;
		vec3_assign_mul(&mass_moment, octant->center.mass);
		return mass_moment;
	}

	mass_moment = zero_vec;
	mass_t total_mass = 0.0;

	for (unsigned i = 0; i < octant_children; i++) {
		const octant_id_t id = octant->children[i];
		struct octant *child = octant_arena_resolve(&octant_arena, id);
		if (unlikely(child == nullptr))
			continue;

		const struct vec3 child_center = octant_update_center(child);
		vec3_assign_add(&mass_moment, &child_center);
		total_mass += child->center.mass;
	}

	octant->center.mass = total_mass;
	octant->center.pos = mass_moment;
	vec3_assign_div(&octant->center.pos, octant->center.mass);

	return mass_moment;
}

extern void
octant_update_accel(const struct octant *octant, const struct point_mass *pm,
	scalar_t len, scalar_t theta_threshold, struct vec3 *accel)
{
	struct vec3 diff = octant->center.pos;
	vec3_assign_sub(&diff, &pm->pos);

	if (octant_is_leaf(octant)) {
		if (unlikely(!vec3_eql(&diff, &zero_vec))) {
			const scalar_t distance = vec3_dist(&diff, &zero_vec);
			vec3_assign_muladd(accel, &diff, distance);
		}

		return;
	}

	const scalar_t distance_sq = vec3_dist_sq(&diff, &zero_vec);
	if (distance_sq > theta_threshold) {
		const scalar_t distance = sqrt(distance_sq);
		vec3_assign_muladd(accel, &diff, distance);

		return;
	}

	for (unsigned i = 0; i < octant_children; i++) {
		if (unlikely(octant->children[i] == invalid_id))
			continue;

		const octant_id_t id = octant->children[i];
		const struct octant *child = octant_arena_resolve(&octant_arena, id);
		octant_update_accel(child, pm, len / 2.0, theta_threshold / 4.0, accel);
	}
}

static void
octant_init(struct octant *octant, const struct point_mass *center,
	struct vec3 corner)
{
	octant->center = *center;
	octant->corner = corner;
	octant->flags = OCTANT_IS_LEAF;
	memset(octant->children, 0, sizeof(octant->children));
}

static void
octant_insert_child(struct octant *octant, const struct point_mass *pm,
	scalar_t child_len)
{
	const struct octant_insertion insertion
		= octant_insertion_prepare(octant, pm, child_len);
	const octant_id_t id = octant->children[insertion.children_idx];
	struct octant *child = octant_arena_resolve(&octant_arena, id);
	if (child != nullptr)
		octant_insert(child, pm, child_len);
	else
		octant_insert_leaf(octant, pm, &insertion);
}

static void
octant_insert_leaf(struct octant *octant, const struct point_mass *pm,
	const struct octant_insertion *insertion)
{
	const octant_id_t id = octant_arena_alloc(&octant_arena);
	struct octant *child = octant_arena_resolve(&octant_arena, id);
	if (unlikely(child == nullptr))
		abort_oom();

	octant->children[insertion->children_idx] = id;
	octant_init(child, pm, insertion->corner);
}

static inline bool
octant_is_leaf(const struct octant *octant)
{
	const enum octant_flags_t flags
		= atomic_load_explicit(&octant->flags, memory_order_relaxed);
	return flags == OCTANT_IS_LEAF;
}

static inline bool
octant_absorbs(const struct octant *octant, const struct point_mass *pm,
	scalar_t child_len)
{
	return scalar_eql(child_len, 0.0) || vec3_eql(&octant->center.pos, &pm->pos);
}

// FIXME: declare first - then put the threadsafe variants here
//
static void
octant_insert_child_threadsafe(struct octant *octant,
	const struct point_mass *pm, scalar_t child_len,
	struct octant_cache *cache)
{
	const struct octant_insertion insertion = octant_insertion_prepare(
		octant, pm, child_len);

	atomic octant_id_t *ptr = (atomic octant_id_t *)&octant->children[insertion.children_idx];
	octant_id_t id = atomic_load_explicit(ptr, memory_order_acquire);
	struct octant *child = octant_arena_resolve(&octant_arena, id);
	if (child == nullptr) {
		id =
	}

	octant_insert_threadsafe(octant, pm, child_len, cache);
}

static octant_id_t
octant_insert_leaf_threadsafe(struct octant *octant,
	const struct point_mass *pm, const struct octant_insertion *insertion,
	struct octant_cache *cache)
{
	const octant_id_t id = octant_cache_alloc(cache);
	struct octant *child = octant_arena_resolve(&octant_arena, id);
	if (unlikely(child == nullptr))
		abort_oom();

	atomic octant_id_t *ptr = (atomic octant_id_t *)(&octant->children[insertion->children_idx]);
	octant_id_t expected = invalid_id;
	const bool success = atomic_compare_exchange_strong_explicit(
		ptr, &expected, id,
		memory_order_release, memory_order_acquire);

	if (likely(success))
		return invalid_id;

	/*
	 * If the CAS failed and another non-null octant ID was found at the
	 * children array position, another thread must have already inserted an
	 * octant there. In this case the unsuccesful thread needs to abort and
	 * instead insert its current particle within that found octant.
	 */
	octant_cache_free_last(cache);
	return expected;
}

static inline bool
octant_try_convert_leaf_threadsafe(struct octant *octant)
{
	enum octant_flags_t flags = atomic_load_explicit(&octant->flags,
		memory_order_relaxed);

	while (true) {
		if ((flags & OCTANT_IS_NODE) != 0)
			return false;

		if (flags == OCTANT_IS_BUSY) {
			cpu_relax();
			flags = atomic_load_explicit(&octant->flags, memory_order_relaxed);
			continue;
		}

		const bool success = atomic_compare_exchange_weak_explicit(
			&octant->flags, &flags, OCTANT_IS_NODE, memory_order_acquire,
			memory_order_relaxed);

		if (likely(success))
			return true;

		cpu_relax();
	}
}

static inline bool
octant_try_absorb_threadsafe(struct octant *octant, const struct point_mass *pm,
	scalar_t child_len)
{
	enum octant_flags_t flags = atomic_load_explicit(&octant->flags,
		memory_order_relaxed);

	if ((flags & OCTANT_IS_NODE) != 0)
		return false;
	if (!octant_absorbs(octant, pm, child_len))
		return false;

	while (true) {
		if (flags == OCTANT_IS_BUSY) {
			cpu_relax();
			flags = atomic_load_explicit(&octant->flags, memory_order_relaxed);
			continue;
		}

		const bool success = atomic_compare_exchange_weak_explicit(
			&octant->flags, &flags, OCTANT_IS_BUSY, memory_order_acquire,
			memory_order_relaxed);

		if (likely(success))
			break;

		if ((flags & OCTANT_IS_NODE) != 0)
			return false;
	}

	octant->center.mass += pm->mass;
	atomic_store_explicit(&octant->flags, OCTANT_IS_LEAF, memory_order_release);

	return true;
}

static struct octant_insertion
octant_insertion_prepare(const struct octant *octant,
	const struct point_mass *pm, scalar_t child_len)
{
	struct octant_insertion insertion;

	const scalar_t xe = octant->corner.x + child_len;
	if (pm->pos.x > xe) {
		insertion.corner.x = xe;
		insertion.children_idx = 1;
	} else {
		insertion.corner.x = octant->corner.x;
		insertion.children_idx = 0;
	}

	const scalar_t ye = octant->corner.y + child_len;
	if (pm->pos.y > ye) {
		insertion.corner.y = ye;
		insertion.children_idx += 2;
	} else
		insertion.corner.y = octant->corner.y;

	const scalar_t ze = octant->corner.z + child_len;
	if (pm->pos.z > ze) {
		insertion.corner.z = ze;
		insertion.children_idx += octant_children / 2;
	} else
		insertion.corner.y = octant->corner.y;

	return insertion;
}

/* Octant arena */

static inline void
octant_arena_reset(struct octant_arena *arena)
{
	for (unsigned i = 0; i < arena->node_count; i++)
		arena->nodes[i].idx = arena->nodes[i].start;
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
#ifdef OPTION_LARGE_SIMULATION
	return id;
#else
	assert(id == invalid_octant || id < ocant_count, "invalid octant ID");
	if (id == invalid_id)
		return (nullptr);
	return &arena->memory[id];
#endif /* OPTION_LARGE_SIMULATION */
}

/* Vector */

static inline struct vec3
vec3_scalar(scalar_t s)
{
	return (struct vec3) { s, s, s };
}

static inline bool
vec3_eql(const struct vec3 *restrict v, const struct vec3 *restrict u)
{
	struct vec3 d = *v;
	vec3_assign_sub(&d, u);

	return vec3_dist_sq(&d, &zero_vec) <= scalar_sq(EPSILON);
}

static inline scalar_t
vec3_dist(const struct vec3 *restrict v, const struct vec3 *restrict u)
{
	return sqrt(vec3_dist_sq(v, u));
}

static inline scalar_t
vec3_dist_sq(const struct vec3 *restrict v, const struct vec3 *restrict u)
{
	assert(v != u, "invalid vector aliasing");
	return scalar_sq(v->x - u->x) + scalar_sq(v->y - u->y)
		+ scalar_sq(v->z - u->z);
}

static inline void
vec3_assign_add(struct vec3 *restrict v, const struct vec3 *restrict u)
{
	assert(v != u, "invalid vector aliasing");
	v->x += u->x;
	v->y += u->y;
	v->z += u->z;
}

static inline void
vec3_assign_sub(struct vec3 *restrict v, const struct vec3 *restrict u)
{
	assert(v != u, "invalid vector aliasing");
	v->x -= u->x;
	v->y -= u->y;
	v->z -= u->z;
}

static inline void
vec3_assign_mul(struct vec3 *v, scalar_t t)
{
	v->x *= t;
	v->y *= t;
	v->z *= t;
}

static inline void
vec3_assign_div(struct vec3 *v, scalar_t t)
{
	assert(t != 0.0, NULL);
	v->x /= t;
	v->y /= t;
	v->z /= t;
}

static inline void
vec3_assign_muladd(struct vec3 *restrict v, const struct vec3 *restrict u,
	scalar_t t)
{
	assert(v != u, "invalid vector aliasing");
	v->x = fma(u->x, t, v->x);
	v->y = fma(u->y, t, v->y);
	v->z = fma(u->z, t, v->z);
}

/* Scalar */

static inline bool
scalar_eql(scalar_t a, scalar_t b)
{
	return (fabs(a - b) <= EPSILON);
}

static inline scalar_t
scalar_sq(scalar_t x)
{

	return (x * x);
}

// FIXME: superfluous
static inline scalar_t
scalar_max(scalar_t a, scalar_t b)
{
	return fmax(a, b);
}

// FIXME: superfluos
static inline scalar_t
scalar_sqrt(scalar_t x)
{
	return sqrt(x);
}

// FIXME: superfluous
static inline scalar_t
scalar_abs(scalar_t x)
{
	return fabs(x);
}
