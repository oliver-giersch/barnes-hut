#define _XOPEN_SOURCE 700

#include "barnes-hut/phys.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

#include <errno.h>
#include <math.h>
#include <string.h>

#include "barnes-hut/arena.h"
#include "barnes-hut/common.h"
#include "barnes-hut/options.h"

// The zero/origin vector.
static const struct vec3 zero_vec = { 0.0, 0.0, 0.0 };

// Returns `x * x`.
static inline float
sq(float x)
{
	return x * x;
}

// Returns `true` if a and b are equal.
static inline bool
feql(float a, float b)
{
	static const float eps = 0.001;
	return fabsf(a - b) <= eps;
}

// Returns a random float between 0.0 and 1.0.
static inline float
randomf(void)
{
	return (float)random() / (float)RAND_MAX;
}

// Returns the morton number for the given x, y, z coordinates.
static inline uint64_t morton_number(unsigned x, unsigned y, unsigned z);
static inline int sort_by_z_curve(const struct accel_particle *p0,
	const struct accel_particle *p1);
static struct vec3 gforce(const struct particle *p0, const struct particle *p1);

// Adds vector `u` to vector `v`.
static inline void vec3_addassign(struct vec3 *v, const struct vec3 *u);
// Subtracts vector `u` from vector `v`.
static inline void vec3_subassign(struct vec3 *v, const struct vec3 *u);
// Multiplies vector `v` with scalar `t`.
static inline void vec3_mulassign(struct vec3 *v, float t);
// Divides vector `v` by scalar `t`.
static inline void vec3_divassign(struct vec3 *v, float t);
// Returns `true` if vectors `v` and `u` are equal.
static inline bool vec3_eql(const struct vec3 *v, const struct vec3 *u);
// Returns the squared distance between vectors `v` and `u`.
static inline float vec3_dist_sq(const struct vec3 *v, const struct vec3 *u);
// Returns the distance between vectors `v` and `u`.
static inline float vec3_dist(const struct vec3 *v, const struct vec3 *u);

void
randomize_particles(struct accel_particle particles[], float r)
{
	for (size_t p = 0; p < options.particles; p++) {
		const float x	 = randomf() * 2 * r - r;
		const float ymax = sqrtf(sq(r) - sq(x));
		const float y	 = randomf() * 2 * ymax - ymax;
		const float zmax = sqrt(sq(r) - sq(x) - sq(y));
		const float z	 = (!options.flat) ? randomf() * 2 * zmax - zmax : 0.0;

		particles[p] = (struct accel_particle) {
			.part = {
				.pos = { x, y, z },
				.mass = options.max_mass,
			},
			.vel = zero_vec,
		};
	}
}

void
sort_particles(struct accel_particle particles[])
{
	typedef int (*cmp_fn)(const void *, const void *);
	qsort(particles, options.particles, sizeof(struct accel_particle),
		(cmp_fn)&sort_by_z_curve);
}

struct octant_malloc_return_t {
	arena_item_t item;
	struct octant *octant;
};

// Returns `true` if the octant represents a leaf in a particle tree.
static inline bool octant_is_leaf(const struct octant *oct);
// Arena-allocates and initializes a new octant for the given center and
// dimensions.
static inline struct octant_malloc_return_t octant_malloc(
	struct particle center, float x, float y, float z, float len);
// Inserts the given particle into one of the octant's children.
static int octant_insert(struct octant *oct, const struct particle *part);
// Inserts the particle into the given child octant.
static int octant_insert_child(struct octant *oct, const struct particle *part);
// Recursively updates the center point of the given octant.
static struct vec3 octant_update_center(struct octant *oct);
// Recursively updates and applies gravitational force to all particles
// contained in the given octant.
static void octant_update_force(const struct octant *oct,
	const struct particle *part, struct vec3 *force);

int
particle_tree_build(struct particle_tree *tree,
	const struct accel_particle particles[], float radius)
{
	int res;

	if (likely(tree->root != ARENA_NULL))
		arena_reset(&arena);

	// Initialize the root octant covering the entire galaxy.
	struct octant_malloc_return_t root = octant_malloc(particles[0].part,
		-1 * radius, -1 * radius, -1 * radius, 2 * radius);
	if ((tree->root = root.item) == ARENA_NULL)
		return ENOMEM;

	// Insert each remaining particle into the tree.
	for (size_t i = 1; i < options.particles; i++) {
		const struct particle *part = &particles[i].part;
		if ((res = unlikely(octant_insert(root.octant, part))))
			return res;
	}

	octant_update_center(root.octant);
	return 0;
}

float
particle_tree_simulate(const struct particle_tree *tree,
	const struct particle_slice *slice)
{
	struct octant *root = arena_get(&arena, tree->root);
	float max_dist_sq	= 0.0;
	float dist_sq		= 0.0;

	for (size_t p = 0; p < slice->len; p++) {
		struct vec3 force		  = zero_vec;
		struct accel_particle *ap = &slice->from[p];
		octant_update_force(root, &ap->part, &force);

		// Apply the calculated force to the particle's velocity.
		vec3_mulassign(&force, options.dt / ap->part.mass);
		vec3_addassign(&ap->vel, &force);
		// Apply the calculated velocity the particle's position.
		struct vec3 vel_dampened = ap->vel;
		vec3_mulassign(&vel_dampened, options.dt);
		vec3_addassign(&ap->part.pos, &vel_dampened);

		dist_sq = vec3_dist_sq(&zero_vec, &ap->part.pos);
		if (dist_sq > max_dist_sq)
			max_dist_sq = dist_sq;
	}

	return sqrtf(max_dist_sq);
}

static inline bool
octant_is_leaf(const struct octant *oct)
{
	return oct->bodies == 1;
}

static inline struct octant_malloc_return_t
octant_malloc(struct particle center, float x, float y, float z, float len)
{
	const arena_item_t item = arena_malloc(&arena, sizeof(struct octant));
	if (unlikely(item == ARENA_NULL))
		return (struct octant_malloc_return_t) { ARENA_NULL, NULL };

	struct octant *oct = arena_get(&arena, item);

	*oct = (struct octant) {
		.center = center,
		.x		= x,
		.y		= y,
		.z		= z,
		.len	= len,
		.bodies = 1,
	};

	for (unsigned c = 0; c < OTREE_CHILDREN; c++)
		oct->children[c] = ARENA_NULL;

	return (struct octant_malloc_return_t) { item, oct };
}

static int
octant_insert(struct octant *oct, const struct particle *part)
{
	int res;

	if (octant_is_leaf(oct)) {
		const bool absorb = vec3_eql(&oct->center.pos, &part->pos)
			|| feql(oct->len / 2.0, 0.0);
		if (absorb) {
			oct->center.mass += part->mass;
			return 0;
		}

		if ((res = octant_insert_child(oct, &oct->center)))
			return res;
	}

	oct->center.mass += part->mass;
	if ((res = octant_insert_child(oct, part)))
		return res;

	return 0;
}

static int
octant_insert_child(struct octant *oct, const struct particle *part)
{
	const float sub_len = oct->len / 2.0;

	float x	   = oct->x;
	float y	   = oct->y;
	float z	   = oct->z;
	unsigned c = 0;

	// Determine, if pos lies in left (0/2) or right (1/3) octant.
	if (part->pos.x > x + sub_len) {
		c = 1;
		x += sub_len;
	}
	// Determine, if pos lies in bottom (0/1) or top (2/3) octant.
	if (part->pos.y > y + sub_len) {
		c += 2;
		y += sub_len;
	}
	// Determine, if pos lies in front or back octant.
	if (part->pos.z > z + sub_len) {
		c += (OTREE_CHILDREN / 2);
		z += sub_len;
	}

	if (oct->children[c] != ARENA_NULL) {
		struct octant *child = arena_get(&arena, oct->children[c]);
		return octant_insert(child, part);
	}

	struct octant_malloc_return_t child
		= octant_malloc(*part, x, y, z, sub_len);
	if (unlikely(child.item == ARENA_NULL))
		return ENOMEM;
	oct->children[c] = child.item;

	return 0;
}

static struct vec3
octant_update_center(struct octant *oct)
{
	struct vec3 new_center;
	if (octant_is_leaf(oct)) {
		new_center = oct->center.pos;
		vec3_mulassign(&new_center, oct->center.mass);
		return new_center;
	}

	new_center = zero_vec;
	for (unsigned c = 0; c < OTREE_CHILDREN; c++) {
		if (oct->children[c] != ARENA_NULL) {
			struct octant *child = arena_get(&arena, oct->children[c]);
			const struct vec3 child_center = octant_update_center(child);
			vec3_addassign(&new_center, &child_center);
		}
	}

	oct->center.pos = new_center;
	vec3_divassign(&oct->center.pos, oct->center.mass);

	return new_center;
}

static void
octant_update_force(const struct octant *oct, const struct particle *part,
	struct vec3 *force)
{
	if (octant_is_leaf(oct)) {
		if (!vec3_eql(&oct->center.pos, &part->pos)) {
			const struct vec3 gf = gforce(part, &oct->center);
			vec3_addassign(force, &gf);
		}

		return;
	}

	const float radius = vec3_dist(&part->pos, &oct->center.pos);
	if (oct->len / radius < options.theta) {
		const struct vec3 gf = gforce(part, &oct->center);
		vec3_addassign(force, &gf);
	} else {
		for (unsigned c = 0; c < OTREE_CHILDREN; c++)
			if (oct->children[c] != ARENA_NULL) {
				const struct octant *child
					= arena_get(&arena, oct->children[c]);
				octant_update_force(child, part, force);
			}
	}
}

static inline void
vec3_addassign(struct vec3 *v, const struct vec3 *u)
{
	v->x += u->x;
	v->y += u->y;
	v->z += u->z;
}

static inline void
vec3_subassign(struct vec3 *v, const struct vec3 *u)
{
	v->x -= u->x;
	v->y -= u->y;
	v->z -= u->z;
}

static inline void
vec3_mulassign(struct vec3 *v, float t)
{
	v->x *= t;
	v->y *= t;
	v->z *= t;
}

static inline void
vec3_divassign(struct vec3 *v, float t)
{
	vec3_mulassign(v, 1 / t);
}

static inline bool
vec3_eql(const struct vec3 *v, const struct vec3 *u)
{
	return feql(v->x, u->x) && feql(v->y, u->y) && feql(v->z, u->z);
}

static inline float
vec3_dist_sq(const struct vec3 *v, const struct vec3 *u)
{
	return sq(v->x - u->x) + sq(v->y - u->y) + sq(v->z - u->z);
}

static inline float
vec3_dist(const struct vec3 *v, const struct vec3 *u)
{
	return sqrtf(vec3_dist_sq(v, u));
}

static inline uint64_t
morton_number(unsigned x, unsigned y, unsigned z)
{
	uint64_t res = 0;

	for (uint64_t i = 0; i < sizeof(uint64_t) * 8 / 3; i++)
		res |= ((x & ((uint64_t)1 << i)) << 2 * i)
			| ((y & ((uint64_t)1 << i)) << (2 * i + 1))
			| ((z & ((uint64_t)1 << i)) << (2 * i + 2));

	return res;
}

static inline int
sort_by_z_curve(const struct accel_particle *p0,
	const struct accel_particle *p1)
{
	const unsigned x0 = (unsigned)p0->part.pos.x;
	const unsigned y0 = (unsigned)p0->part.pos.y;
	const unsigned z0 = (unsigned)p0->part.pos.z;
	const unsigned x1 = (unsigned)p1->part.pos.x;
	const unsigned y1 = (unsigned)p1->part.pos.y;
	const unsigned z1 = (unsigned)p1->part.pos.z;

	const uint64_t m0 = morton_number(x0, y0, z0);
	const uint64_t m1 = morton_number(x1, y1, z1);
	if (m0 < m1)
		return -1;
	else if (m0 > m1)
		return 1;

	return 0;
}

static struct vec3
gforce(const struct particle *p0, const struct particle *p1)
{
	static const float G		= 6.6726e-11;
	static const float min_dist = 2.0;

	if (unlikely(vec3_eql(&p0->pos, &p1->pos)))
		return zero_vec;

	float dist = vec3_dist(&p0->pos, &p1->pos);
	if (dist < min_dist)
		dist = min_dist;

	const float qd = dist * dist * dist;
	const float gm = G * p0->mass * p1->mass;

	struct vec3 result = p1->pos;
	vec3_subassign(&result, &p0->pos);
	vec3_mulassign(&result, gm / qd);

	return result;
}
