#include "../include/phys.h"

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

#include <math.h>

#include "../include/arena.h"
#include "../include/common.h"
#include "../include/options.h"

static inline float
sq(float x)
{
	return x * x;
}

static inline bool
feql(float a, float b)
{
	static const float eps = 0.0001;
	return fabsf(a - b) <= eps;
}

const struct vec2 zvec = { 0.0, 0.0 };
const struct vec2 uvec = { 1.0, 1.0 };

static inline void vec2_addassign(struct vec2 *v, const struct vec2 *u);
static inline void vec2_subassign(struct vec2 *v, const struct vec2 *u);
static inline void vec2_mulassign(struct vec2 *v, float t);
static inline void vec2_divassign(struct vec2 *v, float t);
static inline bool vec2_eql(const struct vec2 *v, const struct vec2 *u);
static inline float vec2_dist_sq(const struct vec2 *v, const struct vec2 *u);
static inline float vec2_dist(const struct vec2 *v, const struct vec2 *u);
static inline bool vec2_is_contained(const struct vec2 *v, float x, float y,
	float len);

struct quadrant {
	struct particle center;
	float x, y, len;
	unsigned bodies;
	struct quadrant *children[4];
};

static inline bool quadrant_is_leaf(const struct quadrant *quad);
static inline struct quadrant *quadrant_malloc(
	struct arena *, struct particle center, float x, float y, float len);
static int quadrant_insert(struct quadrant *quad, const struct particle *part,
	struct arena *arena);
static int quadrant_insert_child(struct quadrant *quad,
	const struct particle *part, struct arena *arena);
static struct vec2 quadrant_update_center(struct quadrant *quad);
static void quadrant_update_force(struct quadrant *quad,
	const struct particle *part, struct vec2 *force);

struct vec2 gforce(const struct particle *p0, const struct particle *p1);

struct particle_tree {
	struct quadrant *root;
	size_t num_particles;
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

	*tree = (struct particle_tree) {
		.root = NULL,
		.num_particles = options.particles,
	};

	return tree;
}

int
particle_tree_build(struct particle_tree *tree, float radius,
	struct arena *arena)
{
	if (likely(tree->root != NULL))
		arena_reset(arena);

	tree->root = quadrant_malloc(arena, tree->particles[0].part,
		-1 * radius, -1 * radius, 2 * radius);
	if (unlikely(tree->root == NULL))
		return -1;

	for (size_t i = 1; i < tree->num_particles; i++) {
		const struct particle *part = &tree->particles[i].part;
		if (quadrant_insert(tree->root, part, arena))
			return -1;
	}

	quadrant_update_center(tree->root);
	return 0;
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
		.center = center,
		.x = x,
		.y = y,
		.len = len,
		.children = { NULL, NULL, NULL, NULL },
	};

	return quad;
}

static int
quadrant_insert(struct quadrant *quad, const struct particle *part,
	struct arena *arena)
{
	int res = 0;

	if (quadrant_is_leaf(quad)) {
		if (vec2_eql(&quad->center.pos, &part->pos)
			|| feql(quad->len / 2, 0.0)
		) {
			quad->center.mass *= part->mass;
			return res;
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
	float x = quad->x;
	float y = quad->y;
	unsigned c;

	if (vec2_is_contained(&part->pos, x, y, len))
		c = 0;
	else if (vec2_is_contained(&part->pos, x + len, y, len))
		c = 1, x = x + len;
	else if (vec2_is_contained(&part->pos, x, y + len, len))
		c = 2, y = y + len;
	else
		c = 3, x = x + len, y = y + len;

	// check if this recursive call is tail-call optimized!
	if (likely(quad->children[c] != NULL))
		return quadrant_insert(quad->children[c], part, arena);

	quad->children[c] = quadrant_malloc(arena, *part, x, y, len);
	if (unlikely(quad->children[c] == NULL))
		return -1;
	return 0;
}


static void
quadrant_update_force(struct quadrant *quad, const struct particle *part,
	struct vec2 *force)
{
	if (unlikely(quadrant_is_leaf(quad))) {
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
		for (unsigned c = 0; c < 4; c++) {
			if (quad->children[c] != NULL)
				quadrant_update_force(quad->children[c], part, force);
		}
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
	return feql(v->x, u->y) && feql(v->y, u->y);
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
