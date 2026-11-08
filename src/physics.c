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
static void octant_insert_child_threadsafe(struct octant *octant,
	const struct point_mass *pm, scalar_t len);
static inline bool octant_is_leaf(const struct octant *octant);
static inline bool octant_is_leaf_threadsafe(const struct octant *octant);

/* Octant insertion */

struct octant_insertion {
	struct vec3 corner;
	uint8_t children_idx;
};

static struct octant_insertion octant_insertion_prepare(
	const struct octant *octant, const struct point_mass *pm,
	scalar_t child_len);

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
	if (octant_is_leaf(octant)) { }

	/* Recursively insert the particle. */
	octant_insert_child(octant, pm, child_len);
}

extern void
octant_insert_threadsafe(struct octant *octant, const struct point_mass *pm,
	scalar_t len, struct octant_cache *cache)
{
	const scalar_t child_len = len / 2.0;
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

		const struct octant *child
			= octant_arena_resolve(&octant_arena, octant->children[i]);
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

// FIXME: declare first

static inline bool
octant_try_absorb_non_blocking(struct octant *octant, struct point_mass *pm,
	scalar_t child_len)
{
	uint32_t is_node
		= atomic_load_explicit(&octant->flags, memory_order_relaxed);
	if (is_node == OCTANT_IS_NODE)
		return (false);

	while (true) {
		if (unlikely(is_node == OCTANT_IS_ABSORBING)) {
			is_node = __atomic_load_n(&octant->flags, memory_order_relaxed);
			__atomic_store();
			continue;
		}

		enum octant_flags_t expected = OCTANT_IS_LEAF;
		const auto res = __atomic_compare_exchange_weak_n(&octant->is_node,
			&expected, OCTANT_IS_ABSORBING, memory_order_acquire,
			memory_order_relaxed);

		if (res)
			break;
	}
}

// FIXME: end

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
