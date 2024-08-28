#ifndef BARNES_HUT_PHYS_H
#define BARNES_HUT_PHYS_H

#include <stddef.h>

#include "barnes-hut/arena.h"

// A 3-dimensional vector.
struct vec3 {
	float x, y, z;
};

// A point-mass particle.
struct point_mass {
	struct vec3 pos;
	float mass;
};

// A moving point-mass particle.
struct particle {
	struct point_mass part;
	struct vec3 vel;
};

// Randomizes the coordinates of the given list of particles.
void randomize_particles(struct particle part[], float r);
// Sorts the the given list of particles by a Z-curve ordering.
void sort_particles(struct particle part[]);

// A consecutive view into the global array of particles.
struct particle_slice {
	// The slice's offset in the source array.
	size_t offset;
	// The slice's length.
	size_t len;
	// The slice's start pointer.
	struct particle *from;
};

// An eight-way partition of a 3-dimensional space containing particles.
#define OTREE_CHILDREN 8
struct octant {
	// The octant's center point mass (cumulative over all contained bodies).
	struct point_mass center;
	// The octant's dimensions (lower-left corner and width).
	float x, y, z, len;
	// The number of all bodies contained within the octant.
	unsigned bodies;
	// The octant's sub-octants (0-3 are (-z)-coords, 4-7 are (+z)-coords).
	//
	//        |---|---|
	//        | 2 | 3 |
	//        |---|---|
	// |---|---|0 | 1 |
	// | 6 | 7 |--|---|
	// |---|---|
	// | 4 | 5 |
	// |---|---|
	arena_item_t children[OTREE_CHILDREN];
};

// A tree of octants containing particles.
struct particle_tree {
	// The particle tree's root octant.
	arena_item_t root;
};

// Recursively constructs the tree structure for the current simulation step.
int particle_tree_build(struct particle_tree *tree,
	const struct particle particles[], float radius);

// Executes the current simulation step by updating all particles encompassed
// by the given slice.
//
// Returns the furthest distance to the center of all updated particles.
float particle_tree_simulate(const struct particle_tree *tree,
	const struct particle_slice *slice);

#endif // BARNES_HUT_PHYS_H
