#include "barnes-hut/physics.h"

#include <stdatomic.>
#include <tgmath.h>

#include "barnes-hut/def.h"

static const struct vec3 zero_vec = { 0.0, 0.0, 0.0 };

static inline bool vec3_eql(const struct vec3 *v, const struct vec3 *u);
static inline scalar_t vec3_dist(const struct vec3 *v, const struct vec3 *u);
static inline scalar_t vec3_dist_sq(const struct vec3 *v, const struct vec3 *u);

static inline void vec3_assign_add(struct vec3 *v, const struct vec3 *u);
static inline void vec3_assign_mul(struct vec3 *v, const struct vec3 *u);
static inline void vec3_assign_div(struct vec3 *v, const struct vec3 *u);
static inline void vec3_assign_addmul(struct vec3 *v, const struct *u,
	scalar_t t);

static inline scalar_t scalar_sq(scalar_t x);
static inline bool scalar_eql(scalar_t a, scalar_t b);
static inline scalar_t scalar_random(void);
static inline scalar_t scalar_max(scalar_t a, scalar_t b);
static inline scalar_t scalar_sqrt(scalar_t x);
static inline scalar_t scalar_abs(scalar_t x);

// FIXME: should be called by thread pool
struct octant *
octant_tree_build(const struct particle particles[], scalar_t radius, scalar_t len)
{
	const auto corner = -radius;

	const auto id = octant_arena_alloc(&arena);
	if (id == invalid_id)
		abort_oom();
	auto *root = octant_arena_get(&arena, id);
	octant_init(root, particles[0].pm, corner, corner, corner);

	for (size_t i = 1; i < particle_count; i++) {
		const auto *pm = &particles[i].pm;
		octant_insert(root, pm, len);
	}
}

// FIXME: declare first
void
octant_insert_non_blocking(struct octant *octant, struct point_mass *pm, scalar_t len, struct octant_cache *cache)
{
	const auto child_len = len / 2.0;


}

static inline bool
octant_try_absorb_non_blocking(struct octant *octant, struct point_mass *pm,
	scalar_t child_len)
{
	uint32_t is_node = atomic_load_explicit(&octant->is_node,
		memory_order_relaxed);
	if (is_node == OCTANT_IS_NODE)
		return (false);

	while (true) {
		if (unlikely(is_node == OCTANT_IS_ABSORBING)) {
			is_node = atomic_load_explicit(&octant->is_node, memory_order_relaxed);
			continue;
		}

		enum octant_flags_t expected = OCTANT_IS_LEAF;
		const auto res = atomic_compare_exchange_weak_explicit(&octant->is_node,
			&expected, OCTANT_IS_ABSORBING, memory_order_acquire, memory_order_relaxed);

		if (res)
			break;
	}
}

// FIXME: end

static inline void
vec3_assign_add(struct vec3 *u, const struct vec3 *v)
{
	v->x += u->x;
	v->y += u->y;
	v->z += u->z;
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
	return (fmax(a, b,));
}

// FIXME: superfluos
static inline scalar_t
scalar_sqrt(scalar_t x)
{
	return (sqrt(x));
}

// FIXME: superfluous
static inline scalar_t
scalar_abs(scalar_t x)
{
	return (fabs(x));
}
