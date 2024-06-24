#define _XOPEN_SOURCE 700

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <assert.h> // TODO: use custom assert
#include <math.h>
#include <pthread.h>
#include <string.h>

#include <sys/mman.h>

#define likely(cond) __builtin_expect(!!(cond), 1)
#define unlikely(cond) __builtin_expect(!!(cond), 0)
#define aligned(align) __attribute__((aligned((align))))

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

// declare & define in arena.h
struct arena {
	void *memory;
	void *curr;
	void *end;
};

// declare in arena.h
int
arena_init(struct arena *arena, size_t size)
{
	arena->memory = mmap(NULL, size, PROT_READ | PROT_WRITE,
		MAP_PRIVATE | MAP_ANON, -1, 0);
	if (unlikely(arena->memory == MAP_FAILED))
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

// define in arena.h (static inline)
void *
arena_malloc(struct arena *arena, size_t size)
{
	void *item = arena->curr;
	if (unlikely(item > arena->end))
		return NULL;
	arena->curr = (void *)((uintptr_t)arena->curr + size);
	return item;
}

// declare & define in phys.h
struct work_slice {
	struct moving_particle *from;
	size_t from_offset;
	size_t end;
};

// main.c
struct thread_state {
	unsigned id;
	struct arena arena;
	struct particle_tree *tree;
	struct work_slice slice;
	float max_dist; // TODO: use as radius R/W, synchronized
};

// main.c
static struct threads {
	unsigned len;
	struct thread_state states[];
} *tls = NULL;

static struct moving_particles *shared_particles;

static inline struct arena *
thread_arena(unsigned id)
{
	assert(tls != NULL && "TLS uninitialized");
	return &tls->states[id].arena;
}

static struct threads *
init_tls2(void)
{
	return NULL;
}

// TODO: mmap arena memory on each thread for NUMA awareness (firsttouch)
static int
init_tls(unsigned threads, size_t arena_sz)
{
	unsigned id;

	tls = malloc(
		sizeof(struct threads) + (sizeof(struct thread_state) * threads));
	if (unlikely(tls == NULL))
		return -1;
	tls->len = threads;
	for (id = 0; id < threads; id++) {
		tls->states[id].id = id;
		// TODO: handle error, unmap
		if (unlikely(arena_init(&tls->states[id].arena, arena_sz)))
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

// declare & define in phys.h, all fns stay in phys.c
struct vec2 {
	float x, y;
};

static const struct vec2 zvec = { 0.0, 0.0 };
static const struct vec2 uvec = { 1.0, 1.0 };

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

static inline struct vec2
vec2_div(struct vec2 v, float t)
{
	vec2_divassign(&v, t);
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

// phys.c
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

// phys.h (declare only!)
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
	if (unlikely(quad == NULL))
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
	if (unlikely(quadrant_is_leaf(quad))) {
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

// declare in phys.h, define in phys.c
struct particle_tree {
	unsigned tid; // TODO: pass as argument to `particle_tree_build`
	struct quadrant *root;
	float theta, dt; // TODO: take from global options
	size_t num_particles;
	struct moving_particle particles[];
};

static struct particle_tree *
particle_tree_init(unsigned tid)
{
	const size_t size = sizeof(struct particle_tree)
		+ (sizeof(struct moving_particle) * options.particles);
	struct particle_tree *tree = malloc(size);
	if (unlikely(tree == NULL))
		return NULL;

	*tree = (struct particle_tree) {
		.tid		   = tid,
		.root		   = NULL,
		.num_particles = options.particles,
	};

	return tree;
}

void
particle_tree_sort(struct particle_tree *tree)
{
	return;
}

// public API, radius can be taken from thread state?
int
particle_tree_build(struct particle_tree *tree, float radius)
{
	struct arena *arena = thread_arena(tree->tid);

	if (likely(tree->root != NULL))
		arena_reset(arena);

	tree->root = quadrant_malloc_init(arena, tree->particles[0].part,
		-1 * radius, -1 * radius, 2 * radius);
	if (unlikely(tree->root == NULL))
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

// define in phys.c, declare in phys.h, struct work_slice must be defined in phys.h!
float
particle_tree_simulate(struct particle_tree *tree,
	const struct work_slice *slice)
{
	struct vec2 force  = {};
	struct vec2 center = {}; // shouldn't center be (r,r)?
	float max_dist	   = 0.0;
	float dist_squared = 0.0;

	for (size_t p = 0; p < slice->end; p++) {
		struct particle *part = &slice->from[p];
		quadrant_update_force(tree->root, part, options.theta, &force);
		const struct vec2 delta_force = vec2_div(force, part->mass);
		vec2_addassign(&slice->from[p].vel, &delta_force);
	}

	// wait for condvar
	// simulate
	// propagate particle positions
	// pass barrier
}

static pthread_barrier_t barrier;

static int
worker_init(unsigned tid)
{
	static const size_t arena_size = (size_t)1 << 20;

	assert(tls != NULL && "TLS not initialized");
	assert(shared_particles != NULL);
	struct thread_state *state = &tls->states[tid];

	state->id = tid;
	if (arena_init(&state->arena, arena_size))
		return -1;
	if ((state->tree = particle_tree_init(tid)) == NULL)
		return -1;

	size_t work_slice	   = options.particles / options.threads;
	const size_t remainder = options.particles % options.threads;
	const size_t start	   = (size_t)tid * work_slice;
	if (tid == options.threads - 1)
		work_slice += remainder;

	memcpy(&state->tree->particles, shared_particles,
		sizeof(struct moving_particle) * options.particles);

	state->slice = (struct work_slice) {
		.from		 = &state->tree->particles[start],
		.end		 = work_slice,
		.from_offset = start,
	};

	return 0;
}

static void
worker_step(struct thread_state *state, unsigned step)
{
	// pthread_barrier_wait(&barrier);
	// run simulation
	// update shared particle array (only slice)
	// pthread_barrier_wait(&barrier);
	// memcpy global particles back into local copy

	(void)pthread_barrier_wait(&barrier);
	// TODO: adjust scale for display?
	if (false) // options.optimize
		if (step % 10 == 0)
			particle_tree_sort(&state->tree);
	// re-load state->radius?
	if (particle_tree_build(&state->tree, 0.0))
		return; // out of memory
	const float radius	= particle_tree_simulate(&state->tree, &state->slice);
	state->tree->radius = radius;
	// TODO: put into function? maybe not, it only requires work slice to be known
	memcpy(&shared_particles[state->slice.from_offset], state->slice.from,
		sizeof(struct moving_particle) * state->slice.end);

	(void)pthread_barrier_wait(&barrier);
	// all threads have copied their loca^l results
	// now copy all shared particles back to our local copy

	// main only: load all farthest - from - center points
	// calculate new universe radius, propagate to all thread states!
}

void *
barnes_hut_worker(void *thread_arg)
{
	const unsigned tid = (unsigned)((uintptr_t)thread_arg);
	if (worker_init(tid))
		// TODO: set global errno?
		return NULL;

	struct thread_state *state = &tls->states[tid];
	for (unsigned step = 0; step < options.steps; step++) {
		worker_step(state, step);
	}

	return NULL;
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

// stays in main.c!
static inline float
randomf(void)
{
	return (float)random() / (float)RAND_MAX;
}

static struct moving_particle *
init_particles(void)
{
	if (options.seed != 0)
		srandom(options.seed);

	struct moving_particle *particles
		= malloc(sizeof(struct moving_particle) * options.particles);
	if (particles == NULL)
		return NULL;

	// TODO: Universe center is at (r, r)
	// - make all initialization around (0, 0), then translate by adding (r, r)
	// - get random x first, between [-r, r]
	// - for y to be located within the circle, the random y-pos must be between
	//   [-sqrt(r^2 - x^2), sqrt(r^2 - x^2)]
	const float r = options.radius;
	const float d = r * 2;

	for (size_t p = 0; p < options.particles; p++) {
		const float x		   = (randomf() * d) - r;
		const float tmp		   = sqrtf(r * r - x * x);
		const float y		   = (randomf() * 2 * tmp) - tmp;
		particles[p].part.pos  = (struct vec2) { x + r, y + r };
		particles[p].part.mass = options.max_mass;
		particles[p].vel	   = zvec;
	}

	return particles;
}

int
main(int argc, char *argv[argc])
{
	// parse arguments
	// pin pthread_current() to cpu 0
	// initialize random particle list

	shared_particles = init_particles();
	if (shared_particles == NULL)
		return -1;

	tls = init_tls2();
	if (tls == NULL)
		return -1;

	pthread_barrierattr_t attr;
	pthread_barrierattr_init(&attr);
	pthread_barrier_init(&barrier, &attr, options.threads);

	// spawn p - 1 threads
	pthread_t *threads = malloc(sizeof(pthread_t) * options.threads);
	if (threads == NULL)
		return -1;

	unsigned t;
	for (t = 0; t < options.threads + 1; t++) {
		pthread_attr_t attrs;
		pthread_attr_init(&attrs);
		if (pthread_create(&threads[t], &attrs, barnes_hut_worker,
				(void *)(t + 1)))
			goto kill;
		pthread_attr_destroy(&attrs);
	}

	if (worker_init(0))
		goto kill;

	for (unsigned step = 0; step < options.steps; step++) {
		worker_step(0, step);
	}

kill:
	for (unsigned i = 0; i < t; i++)
		pthread_kill(threads[i], SIGKILL);

	return 0;
}
