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
	static const float eps = 0.0001;
	return fabsf(a - b) <= eps;
}

// Returns a random float between 0.0 and 1.0.
static inline float
randomf(void)
{
	return (float)random() / (float)RAND_MAX;
}

static inline uint64_t morton_number(unsigned x, unsigned y, unsigned z);
static inline int sort_by_z_curve(const struct accel_particle *p0,
	const struct accel_particle *p1);
static struct vec3 gforce(const struct particle *p0, const struct particle *p1);

static inline void vec3_addassign(struct vec3 *v, const struct vec3 *u);
static inline void vec3_subassign(struct vec3 *v, const struct vec3 *u);
static inline void vec3_mulassign(struct vec3 *v, float t);
static inline void vec3_divassign(struct vec3 *v, float t);
static inline bool vec3_eql(const struct vec3 *v, const struct vec3 *u);
static inline float vec3_dist_sq(const struct vec3 *v, const struct vec3 *u);
static inline float vec3_dist(const struct vec3 *v, const struct vec3 *u);

void
randomize_particles(struct accel_particle particles[], float r)
{
	for (size_t p = 0; p < options.particles; p++) {
		const float x	 = randomf() * 2 * r - r;
		const float ymax = sqrtf(sq(r) - sq(x));
		const float y	 = randomf() * 2 * ymax - ymax;
		const float zmax = sqrt(sq(r) - sq(x) - sq(y));
		const float z	 = randomf() * 2 * zmax - zmax;

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

#define OTREE_CHILDREN 8
struct octant {
	// The octant's center point mass (cumulative over all contained bodies).
	struct particle center;
	// The octant's dimensions (lower-left corner and width).
	float x, y, z, len;
	// The number of all bodies contained within the octant.
	unsigned bodies;
	// The octant's sub-octants (0-3 are (-z)-coords, 4-7 are (+z)-coords).
	//
	// |---|---|
	// | 2 | 3 |
	// |---|---|
	// | 0 | 1 |
	// |---|---|
	struct octant *children[OTREE_CHILDREN];
};

// Returns `true` if the octant represents a leaf in a particle tree.
static inline bool octant_is_leaf(const struct octant *oct);
// Arena-allocates and initializes a new octant for the given center and
// dimensions.
static inline struct octant *octant_malloc(struct arena *,
	struct particle center, float x, float y, float z, float len);
// Inserts the given particle into one of the octant's children.
static int octant_insert(struct octant *oct, const struct particle *part,
	struct arena *arena);
// Inserts the particle into the given child octant.
static int octant_insert_child(struct octant *oct, const struct particle *part,
	struct arena *arena);
// Recursively updates the center point of the given octant.
static struct vec3 octant_update_center(struct octant *oct);
// Recursively updates and applies gravitational force to all particles
// contained in the given octant.
static void octant_update_force(const struct octant *oct,
	const struct particle *part, struct vec3 *force);

int
particle_tree_build(struct particle_tree *tree,
	const struct accel_particle particles[], float radius, struct arena *arena)
{
	int res;

	if (likely(tree->root != NULL))
		arena_reset(arena);

	// Initialize the root octant covering the entire galaxy.
	tree->root = octant_malloc(arena, particles[0].part, -1 * radius,
		-1 * radius, -1 * radius, 2 * radius);
	if (unlikely(tree->root == NULL))
		return ENOMEM;

	// Insert each remaining particle into the tree.
	for (size_t i = 1; i < options.particles; i++) {
		const struct particle *part = &particles[i].part;
		if ((res = unlikely(octant_insert(tree->root, part, arena))))
			return res;
	}

	octant_update_center(tree->root);
	return 0;
}

float
particle_tree_simulate(const struct particle_tree *tree,
	const struct particle_slice *slice)
{
	float max_dist	   = 0.0;
	float dist_squared = 0.0;

	for (size_t p = 0; p < slice->len; p++) {
		struct vec3 force		  = zero_vec;
		struct accel_particle *ap = &slice->from[p];
		octant_update_force(tree->root, &ap->part, &force);

		// Apply the calculated force to the particle's velocity.
		vec3_divassign(&force, ap->part.mass);
		vec3_addassign(&ap->vel, &force);
		// Apply the calculated velocity the particle's position.
		vec3_addassign(&ap->part.pos, &ap->vel);

		dist_squared = vec3_dist_sq(&zero_vec, &ap->part.pos);
		if (dist_squared > max_dist)
			max_dist = dist_squared;
	}

	return sqrtf(max_dist);
}

static inline bool
octant_is_leaf(const struct octant *oct)
{
	return oct->bodies == 1;
}

static inline struct octant *
octant_malloc(struct arena *arena, struct particle center, float x, float y,
	float z, float len)
{
	struct octant *oct = arena_malloc(arena, sizeof(struct octant));
	if (unlikely(oct == NULL))
		return NULL;

	*oct = (struct octant) {
		.center	  = center,
		.x		  = x,
		.y		  = y,
		.z		  = z,
		.len	  = len,
		.bodies	  = 1,
		.children = { NULL },
	};

	return oct;
}

static int
octant_insert(struct octant *oct, const struct particle *part,
	struct arena *arena)
{
	int res;

	if (octant_is_leaf(oct)) {
		const bool absorb
			= vec3_eql(&oct->center.pos, &part->pos) || feql(oct->len / 2, 0.0);
		if (absorb) {
			oct->center.mass += part->mass;
			return 0;
		}

		if ((res = octant_insert_child(oct, &oct->center, arena)))
			return res;
	}

	oct->center.mass += part->mass;
	if ((res = octant_insert_child(oct, part, arena)))
		return res;

	return 0;
}

static int
octant_insert_child(struct octant *oct, const struct particle *part,
	struct arena *arena)
{
	const float len = oct->len / 2.0;
	float x			= oct->x;
	float y			= oct->y;
	float z			= oct->z;
	unsigned c;

	// Determine, if pos lies in left (0/2) or right (1/3) octant.
	if (part->pos.x <= x + len)
		c = 0;
	else
		c = 1, x += len;
	// Determine, if pos lies in bottom (0/1) or top (2/3) octant.
	if (part->pos.y > y + len)
		c += 2, y += len;
	// Determine, if pos lies in front or back octant.
	if (part->pos.z > z + len)
		c += (OTREE_CHILDREN / 2), z += len;

	if (oct->children[c] != NULL)
		return octant_insert(oct->children[c], part, arena);

	oct->children[c] = octant_malloc(arena, *part, x, y, z, len);
	if (unlikely(oct->children[c] == NULL))
		return ENOMEM;

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
		if (oct->children[c] != NULL) {
			const struct vec3 child_center
				= octant_update_center(oct->children[c]);
			vec3_addassign(&new_center, &child_center);
		}
	}

	vec3_divassign(&new_center, oct->center.mass);
	oct->center.pos = new_center;
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
			if (oct->children[c] != NULL)
				octant_update_force(oct->children[c], part, force);
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

	const float rp3 = dist * dist * dist;
	const float gm	= G * p0->mass * p1->mass;

	struct vec3 result = p1->pos;
	vec3_subassign(&result, &p0->pos);
	vec3_mulassign(&result, gm / rp3);

	return result;
}
