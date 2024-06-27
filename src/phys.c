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
static const struct vec2 zero_vec = { 0.0, 0.0 };

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

static inline void vec2_addassign(struct vec2 *v, const struct vec2 *u);
static inline void vec2_subassign(struct vec2 *v, const struct vec2 *u);
static inline void vec2_mulassign(struct vec2 *v, float t);
static inline void vec2_divassign(struct vec2 *v, float t);
static inline bool vec2_eql(const struct vec2 *v, const struct vec2 *u);
static inline float vec2_dist_sq(const struct vec2 *v, const struct vec2 *u);
static inline float vec2_dist(const struct vec2 *v, const struct vec2 *u);

void
moving_particle_randomize(struct moving_particle *part, float r)
{
	const float x	 = randomf() * 2 * r - r;
	const float ymax = sqrtf(sq(r) - sq(x));
	const float y	 = randomf() * 2 * ymax - ymax;

	*part = (struct moving_particle) {
		.part = {
			.pos = { x, y },
			.mass = options.max_mass,
		},
		.vel = zero_vec,
	};
}

#define QTREE_CHILDREN 4
struct quadrant {
	// The quadrant's center point mass (cumulative over all contained bodies).
	struct particle center;
	// The quadrant's dimensions (lower-left corner and width).
	float x, y, len;
	// The number of all bodies contained within the quadrant.
	unsigned bodies;
	// The quadrant's sub-quadrants
	//   0: lower-left
	//   1: lower-right
	//   2: upper-left
	//   3: upper-right
	struct quadrant *children[QTREE_CHILDREN];
};

static inline bool quadrant_is_leaf(const struct quadrant *quad);
static inline struct quadrant *quadrant_malloc(struct arena *,
	struct particle center, float x, float y, float len);
static int quadrant_insert(struct quadrant *quad, const struct particle *part,
	struct arena *arena);
static int quadrant_insert_child(struct quadrant *quad,
	const struct particle *part, struct arena *arena);
static struct vec2 quadrant_update_center(struct quadrant *quad);
static void quadrant_update_force(struct quadrant *quad,
	const struct particle *part, struct vec2 *force);

static inline uint64_t morton_number(unsigned x, unsigned y, unsigned z);
static inline int sort_by_z_curve(const struct moving_particle *p0,
	const struct moving_particle *p1);
static struct vec2 gforce(const struct particle *p0, const struct particle *p1);

struct particle_tree {
	// The particle tree's root quadrant.
	struct quadrant *root;
	struct moving_particle particles[];
};

struct particle_tree *
particle_tree_init(void)
{
	const size_t size = sizeof(struct particle_tree)
		+ (sizeof(struct moving_particle) * options.particles);
	struct particle_tree *tree = malloc(size);
	if (unlikely(tree == NULL))
		return NULL;

	*tree = (struct particle_tree) { .root = NULL };

	return tree;
}

void
particle_tree_sort(struct particle_tree *tree)
{
	typedef int (*cmp_fn)(const void *, const void *);
	qsort(tree->particles, options.particles, sizeof(struct moving_particle),
		(cmp_fn)&sort_by_z_curve);
}

int
particle_tree_build(struct particle_tree *tree, float radius,
	struct arena *arena)
{
	int res;

	if (likely(tree->root != NULL))
		arena_reset(arena);

	tree->root = quadrant_malloc(arena, tree->particles[0].part, -1 * radius,
		-1 * radius, 2 * radius);
	if (unlikely(tree->root == NULL))
		return ENOMEM;

	for (size_t i = 1; i < options.particles; i++) {
		const struct particle *part = &tree->particles[i].part;
		if ((res = unlikely(quadrant_insert(tree->root, part, arena))))
			return res;
	}

	quadrant_update_center(tree->root);
	return 0;
}

float
particle_tree_simulate(struct particle_tree *tree,
	const struct particle_slice *slice)
{
	const struct vec2 center = zero_vec;

	struct vec2 force  = zero_vec;
	float max_dist	   = 0.0;
	float dist_squared = 0.0;

	for (size_t p = 0; p < slice->len; p++) {
		struct particle *part = &slice->from[p].part;
		quadrant_update_force(tree->root, part, &force);

		struct vec2 delta_force = force;
		vec2_divassign(&delta_force, part->mass);
		vec2_addassign(&slice->from[p].vel, &delta_force);

		dist_squared = vec2_dist_sq(&center, &part->pos);
		if (dist_squared > max_dist)
			max_dist = dist_squared;

		force = zero_vec;
	}

	return sqrtf(max_dist);
}

struct moving_particle *
particle_tree_particles(struct particle_tree *tree)
{
	return &tree->particles[0];
}

static inline bool
quadrant_is_leaf(const struct quadrant *quad)
{
	return quad->bodies == 1;
}

static inline struct quadrant *
quadrant_malloc(struct arena *arena, struct particle center, float x, float y,
	float len)
{
	struct quadrant *quad = arena_malloc(arena, sizeof(struct quadrant));
	if (unlikely(quad == NULL))
		return NULL;

	*quad = (struct quadrant) {
		.center	  = center,
		.x		  = x,
		.y		  = y,
		.len	  = len,
		.bodies	  = 1,
		.children = { NULL, NULL, NULL, NULL },
	};

	return quad;
}

static int
quadrant_insert(struct quadrant *quad, const struct particle *part,
	struct arena *arena)
{
	int res;

	if (quadrant_is_leaf(quad)) {
		if (vec2_eql(&quad->center.pos, &part->pos)
			|| feql(quad->len / 2, 0.0)) {
			quad->center.mass += part->mass;
			return 0;
		}

		if ((res = quadrant_insert_child(quad, &quad->center, arena)))
			return res;
	}

	quad->center.mass += part->mass;
	if ((res = quadrant_insert_child(quad, part, arena)))
		return res;

	return 0;
}

static int
quadrant_insert_child(struct quadrant *quad, const struct particle *part,
	struct arena *arena)
{
	const float len = quad->len / 2.0;
	float x			= quad->x;
	float y			= quad->y;
	unsigned c;

	// Determine, if pos lies in left (0/2) or right (1/3) quadrant.
	if (part->pos.x <= x + len)
		c = 0;
	else
		c = 1, x += len;
	// Determine, if pos lies in bottom (0/1) or top (2/3) quadrant.
	if (part->pos.y > y + len)
		c += 2, y += len;

	if (quad->children[c] != NULL)
		return quadrant_insert(quad->children[c], part, arena);

	quad->children[c] = quadrant_malloc(arena, *part, x, y, len);
	if (unlikely(quad->children[c] == NULL))
		return ENOMEM;

	return 0;
}

static struct vec2
quadrant_update_center(struct quadrant *quad)
{
	struct vec2 new_center;
	if (quadrant_is_leaf(quad)) {
		new_center = quad->center.pos;
		vec2_mulassign(&new_center, quad->center.mass);
		return new_center;
	}

	new_center = zero_vec;
	for (unsigned c = 0; c < QTREE_CHILDREN; c++) {
		if (quad->children[c] != NULL) {
			const struct vec2 child_center
				= quadrant_update_center(quad->children[c]);
			vec2_addassign(&new_center, &child_center);
		}
	}

	vec2_divassign(&new_center, quad->center.mass);
	quad->center.pos = new_center;
	return new_center;
}

static void
quadrant_update_force(struct quadrant *quad, const struct particle *part,
	struct vec2 *force)
{
	if (quadrant_is_leaf(quad)) {
		if (!vec2_eql(&quad->center.pos, &part->pos)) {
			const struct vec2 gf = gforce(part, &quad->center);
			vec2_addassign(force, &gf);
		}

		return;
	}

	const float radius = vec2_dist(&part->pos, &quad->center.pos);
	if (quad->len / radius < options.theta) {
		const struct vec2 gf = gforce(part, &quad->center);
		vec2_addassign(force, &gf);
	} else {
		for (unsigned c = 0; c < QTREE_CHILDREN; c++)
			if (quad->children[c] != NULL)
				quadrant_update_force(quad->children[c], part, force);
	}
}

static inline void
vec2_addassign(struct vec2 *v, const struct vec2 *u)
{
	v->x += u->x;
	v->y += u->y;
}

static inline void
vec2_subassign(struct vec2 *v, const struct vec2 *u)
{
	v->x -= u->x;
	v->y -= u->y;
}

static inline void
vec2_mulassign(struct vec2 *v, float t)
{
	v->x *= t;
	v->y *= t;
}

static inline void
vec2_divassign(struct vec2 *v, float t)
{
	vec2_mulassign(v, 1 / t);
}

static inline bool
vec2_eql(const struct vec2 *v, const struct vec2 *u)
{
	return feql(v->x, u->x) && feql(v->y, u->y);
}

static inline float
vec2_dist_sq(const struct vec2 *v, const struct vec2 *u)
{
	return sq(v->x - u->x) + sq(v->y - u->y);
}

static inline float
vec2_dist(const struct vec2 *v, const struct vec2 *u)
{
	return sqrtf(vec2_dist_sq(v, u));
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
sort_by_z_curve(const struct moving_particle *p0,
	const struct moving_particle *p1)
{
	const unsigned x0 = (unsigned)p0->part.pos.x;
	const unsigned y0 = (unsigned)p0->part.pos.y;
	const unsigned x1 = (unsigned)p1->part.pos.x;
	const unsigned y1 = (unsigned)p1->part.pos.y;

	return (morton_number(x0, y0, 0) < morton_number(x1, y1, 0)) ? -1 : 1;
}

static struct vec2
gforce(const struct particle *p0, const struct particle *p1)
{
	static const float G		= 6.6726e-11;
	static const float min_dist = 4.0;

	if (unlikely(vec2_eql(&p0->pos, &p1->pos)))
		return zero_vec;

	float dist = vec2_dist(&p0->pos, &p1->pos);
	if (dist < min_dist)
		dist = min_dist;

	const float rp3 = dist * dist * dist;
	const float gm	= G * p0->mass * p1->mass;

	struct vec2 result = p1->pos;
	vec2_subassign(&result, &p0->pos);
	vec2_mulassign(&result, gm / rp3);

	return result;
}
