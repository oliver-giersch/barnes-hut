#ifndef BARNES_HUT_PHYS_H
#define BARNES_HUT_PHYS_H

#include <stddef.h>

// A 2-dimensional vector.
struct vec2 {
	float x, y;
};

// A point-mass particle.
struct particle {
	struct vec2 pos;
	float mass;
};

// A moving point-mass particle.
struct moving_particle {
	struct particle part;
	struct vec2 vel;
};

// Randomizes the coordinates of the given particle within a cirle or radius r
void moving_particle_randomize(struct moving_particle *part, float r);

// A consecutive view into the global array of particles.
struct particle_slice {
	// The slice's offset in the source array.
	size_t offset;
	// The slice's length.
	size_t len;
	// The slice's start pointer.
	struct moving_particle *from;
};

// A four-way partition of a 2-dimensional space containing particles.
struct quadrant;

// An fixed size arena allocator for quadrants.
struct arena;

// A tree of quadrants containing particles.
struct particle_tree;

// Allocates and initializes a particle tree.
//
// Returns `NULL` on failure.
struct particle_tree *particle_tree_init(void);

// Sorts the particle tree's local array of particles.
void particle_tree_sort(struct particle_tree *tree);

// Recursively constructs the tree structure for the current simulation step.
int particle_tree_build(struct particle_tree *tree, float radius,
	struct arena *arena);

// Executes the current simulation step by updating all particles encompassed
// by the given slice.
//
// Returns the furthest distance to the center of all updated particles.
float particle_tree_simulate(struct particle_tree *tree,
	const struct particle_slice *slice);

// Returns the particle tree's local copy of the global particle array.
struct moving_particle *particle_tree_particles(struct particle_tree *tree);

#endif // BARNES_HUT_PHYS_H
