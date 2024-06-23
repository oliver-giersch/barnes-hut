#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <assert.h> // TODO: use custom assert
#include <math.h>
#include <pthread.h>

#include <sys/mman.h>

#define likely(cond) (cond)
#define unlikely(cond) (cond)

static inline float
square(float x)
{
	return x * x;
}

static inline bool
feql(float a, float b)
{
	static const float eps = 0.00001;
	return fabsf(a - b) <= eps;
}

struct arena {
	void *memory;
	void *curr;
	void *end;
};

int
arena_init(struct arena *arena, size_t size)
{
	arena->memory = mmap(NULL, size, PROT_READ | PROT_WRITE,
		MAP_PRIVATE | MAP_ANON, -1, 0);
	if (arena->memory == MAP_FAILED)
		return -1;

	arena->curr = arena->memory;
	arena->end	= (void *)((uintptr_t)arena->memory + size);

	return 0;
}

void
arena_deinit(struct arena *arena)
{
	const size_t size = (uintptr_t)arena->end - (uintptr_t)arena->curr;
	munmap(arena->memory, size);
}

void
arena_reset(struct arena *arena)
{
	arena->curr = arena->memory;
}

void *
arena_malloc(struct arena *arena, size_t size)
{
	void *item = arena->curr;
	if (item > arena->end)
		return NULL;
	arena->curr = (void *)((uintptr_t)arena->curr + size);
	return item;
}

struct thread_state {
	unsigned id;
	struct arena arena;
	struct particle_tree *tree;
};

static struct threads {
	unsigned len;
	struct thread_state states[];
} *tls = NULL;

static inline struct arena *
thread_arena(unsigned id)
{
	assert(tls != NULL && "TLS uninitialized");
	return &tls->states[id].arena;
}

// TODO: mmap arena memory on each thread for NUMA awareness (firsttouch)
static int
init_tls(unsigned threads, size_t arena_sz)
{
	unsigned id;

	tls = malloc(
		sizeof(struct threads) + (sizeof(struct thread_state) * threads));
	if (tls == NULL)
		return -1;
	tls->len = threads;
	for (id = 0; id < threads; id++) {
		tls->states[id].id = id;
		// TODO: handle error, unmap
		if (arena_init(&tls->states[id].arena, arena_sz))
			goto error;
	}

	return 0;

error:
	for (unsigned i = 0; i < id; i++) {
		arena_deinit(&tls->states[i].arena);
	}

	free(tls);
	return -1;
};

struct vec2 {
	float x, y;
};

static inline void
vec2_addassign(struct vec2 *v, const struct vec2 *u)
{
	v->x += u->x;
	v->y += u->y;
}

static inline void
vec2_addassign(struct vec2 *v, const struct vec2 *u)
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

static inline struct vec2
vec2_sub(struct vec2 v, const struct vec2 *u)
{
	vec2_subassign(&v, u);
	return v;
}

static inline struct vec2
vec2_mul(struct vec2 v, float t)
{
	vec2_mulassign(&v, t);
	return v;
}

static inline bool
vec2_eql(const struct vec2 *v, const struct vec2 *u)
{
	return feql(v->x, u->x) && feql(v->y, v->y);
}

static inline float
vec2_dist(const struct vec2 *v, const struct vec2 *u)
{
	return sqrtf(square(v->x - u->x) + square(v->y - u->y));
}

static inline bool
vec2_contained_within(const struct vec2 *v, float x, float y, float len)
{
	return (v->x >= x && v->x <= x + len) //  x
		&& (v->y >= y && v->y <= y + len); // y
}

struct particle {
	struct vec2 pos;
	float mass;
};

struct moving_particle {
	struct particle part;
	struct vec2 vel;
};

static struct vec2
gforce(const struct particle *p0, const struct particle *p1)
{
	static const float G		= 6.6726e-11;
	static const float min_dist = 4.0;

	if (vec2_eql(&p0, &p1))
		return (struct vec2) {};

	float dist = vec2_dist(&p0, &p1);
	if (dist < min_dist)
		dist = min_dist;

	const float rp3 = dist * dist * dist;
	const float gm	= G * p0->mass * p1->mass;

	struct vec2 result = vec2_sub(p1->pos, &p0->pos);
	vec2_mulassign(&result, gm / rp3);

	return result;
}

struct quadrant {
	// The quadrant's mass center.
	struct particle center;
	// The quadrant's dimensions starting from the lower left corner.
	float x, y, len;
	// The quadrant's number of contained bodies.
	unsigned bodies;
	// The quadrant's child quadrants.
	struct quadrant *children[4];
};

static inline bool
quadrant_is_leaf(const struct quadrant *quad)
{
	return quad->bodies == 1;
}

static inline struct quadrant *
quadrant_malloc_init(struct arena *arena, struct particle center, float x,
	float y, float len)
{
	struct quadrant *quad = arena_malloc(arena, sizeof(struct quadrant));
	if (quad == NULL)
		return NULL;

	*quad = (struct quadrant) {
		.center	  = center,
		.x		  = x,
		.y		  = y,
		.len	  = len,
		.children = { NULL, NULL, NULL, NULL },
	};

	return quad;
}

static int quadrant_insert(struct quadrant *quad, const struct particle *part,
	struct arena *arena);

static int
quadrant_insert_child(struct quadrant *quad, const struct particle *part,
	struct arena *arena)
{
	const float x	= quad->x;
	const float y	= quad->y;
	const float len = quad->len / 2.0;

	struct vec2 corner;
	unsigned child;

	if (vec2_contained_within(&part->pos, x, y, len))
		child = 0, corner = (struct vec2) { x, y };
	else if (vec2_contained_within(&part->pos, x + len, y, len))
		child = 1, corner = (struct vec2) { x + len, y };
	else if (vec2_contained_within(&part->pos, x, y + len, len))
		child = 2, corner = (struct vec2) { x, y + len };
	else
		child = 3, corner = (struct vec2) { x + len, y + len };

	if (likely(quad->children[child] != NULL))
		return quadrant_insert(quad->children[child], part, arena);

	quad->children[child]
		= quadrant_malloc_init(arena, *part, corner.x, corner.y, len);
	if (unlikely(quad->children[child] == NULL))
		return -1;
	return 0;
}

static int
quadrant_insert(struct quadrant *quad, const struct particle *part,
	struct arena *arena)
{
	if (quad->bodies == 1) {
		if (vec2_eql(&quad->center, &part->pos) || feql(quad->len / 2, 0.0)) {
			quad->center.mass += part->mass;
			return 0;
		}

		if (quadrant_insert_child(quad, &quad->center, arena))
			return -1;
	}

	quad->center.mass += part->mass;
	if (quadrant_insert_child(quad, part, arena))
		return -1;

	return 0;
}

static struct vec2
quadrant_update_center(struct quadrant *quad)
{
	if (unlikely(quad->bodies <= 1))
		return vec2_mul(quad->center.pos, quad->center.mass);

	struct vec2 cpos = {};
	for (unsigned c = 0; c < 4; c++) {
		if (quad->children[c] != NULL) {
			const struct vec2 ccenter
				= quadrant_update_center(quad->children[c]);
			vec2_addassign(&cpos, &ccenter);
		}
	}

	vec2_divassign(&cpos, quad->center.mass);
	quad->center.pos = cpos;
	return cpos;
}

static struct vec2
quadrant_update_force(struct quadrant *quad, const struct particle *part,
	float theta, struct vec2 *force)
{
	if (quadrant_is_leaf(quad)) {
		if (!vec2_eql(&quad->center, &part->pos)) {
			const struct vec2 gf = gforce(part, &quad->center);
			vec2_addassign(force, &gf);
		}

		return;
	}

	const float radius = vec2_dist(&part->pos, &quad->center);
	if (quad->len / radius < theta) {
		const struct vec2 gf = gforce(part, &quad->center);
		vec2_addassign(force, &gf);
	} else {
		for (unsigned c = 0; c < 4; c++) {
			if (quad->children[c] != NULL)
				quadrant_update_force(quad->children[c], part, theta, force);
		}
	}
}

struct particle_tree {
	unsigned tid; // TODO: pass as argument to `particle_tree_build`
	struct quadrant *root;
	float theta, dt, radius; // TODO: take from global options
	size_t num_particles;
	struct moving_particle particles[];
};

static int
particle_tree_init(struct particle_tree *tree, unsigned tid,
	size_t num_particles)
{
	struct moving_particle *particles
		= malloc(sizeof(struct moving_particle) * num_particles);
	if (particles == NULL)
		return -1;

	*tree = (struct particle_tree) {
		.tid		   = tid,
		.root		   = NULL,
		.num_particles = num_particles,
		.particles	   = particles,
	};

	return 0;
}

int
particle_tree_build(struct particle_tree *tree, float radius)
{
	struct arena *arena = thread_arena(tree->tid);

	if (tree->root != NULL)
		arena_reset(arena);

	tree->root = quadrant_malloc_init(arena, tree->particles[0].part,
		-1 * radius, -1 * radius, 2 * radius);
	if (tree->root == NULL)
		return -1;

	for (size_t i = 0; i < tree->num_particles; i++) {
		const struct particle *part = &tree->particles[i];
		if (quadrant_insert(tree->root, part, arena))
			return -1;
	}

	quadrant_update_center(tree->root);
}

void
particle_tree_propagate(struct particle_tree *tree)
{
	// write positions slice into shared particle area
	// pthread_barrier
}

void
particle_tree_simulate(struct particle_tree *tree)
{
	// wait for condvar
	// simulate
	// propagate particle positions
	// pass barrier
}

static pthread_barrier_t barrier;

static int
worker_step(unsigned tid, unsigned step)
{
	// pthread_barrier_wait(&barrier);
	// run simulation
	// update shared particle array (only slice)
	// pthread_barrier_wait(&barrier);
	// memcpy global particles back into local copy
}

void *
barnes_hut_worker(void *thread_arg)
{
	const unsigned tid = (unsigned)((uintptr_t)thread_arg);
	// TODO: init thread state
	// initialize particle tree, work chunk, etc
}

struct options {
	unsigned steps;
	size_t particles;
	float max_mass;
	float radius;
	float theta;
	unsigned threads;
	unsigned seed;
} options;

static void
print_usage(const char *bin)
{
	fprintf(stderr,
		"Usage: %s [options]\n"
		"Options:\n"
		"-t [STEPS], --time=[STEPS]\n"
		"-n [PARTICLES], --num=[PARTICLES]\n"
		"-m [MASS], --mass=[MASS]\n"
		"-r [RADIUS], --radius=[RADIUS]\n"
		"-p [THREADS], --threads=[THREADS]\n"
		"-s [SEED], --seed=[SEED]\n"
		"-h, --help\n"
#ifdef DISPLAY
		"-d [DIM], --dimension=[DIM]\n"
#endif // DISPLAY
		"--theta\n");
}

int
main(int argc, char *argv[argc])
{
	// parse arguments
	// pin pthread_current() to cpu 0
	// initialize random particle list

	// spawn p - 1 threads
	pthread_t *threads = malloc(sizeof(pthread_t) * options.threads);
	if (threads == NULL)
		return -1;

	unsigned t;
	for (t = 0; t < options.threads; t++) {
		pthread_attr_t attrs;
		pthread_attr_init(&attrs);
		if (pthread_create(&threads[t], &attrs, barnes_hut_worker, (void *)t))
			goto kill;
		pthread_attr_destroy(&attrs);
	}

	for (unsigned step = 0; step < options.steps; step++) {
		// run thread step fn
		//
	}

kill:
	for (unsigned i = 0; i < t; i++)
		pthread_kill(threads[i], SIGKILL);

	return 0;
}
