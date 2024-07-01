#ifndef BARNES_HUT_PHYS_H
#define BARNES_HUT_PHYS_H

#include <stddef.h>

// A 3-dimensional vector.
struct vec3 {
	float x, y, z;
};

// A point-mass particle.
struct particle {
	struct vec3 pos;
	float mass;
};

// A moving point-mass particle.
struct accel_particle {
	struct particle part;
	struct vec3 vel;
};

// Randomizes the coordinates of the given list of particles.
void randomize_particles(struct accel_particle part[], float r);
// Sorts the the given list of particles by a Z-curve ordering.
void sort_particles(struct accel_particle part[]);

// A consecutive view into the global array of particles.
struct particle_slice {
	// The slice's offset in the source array.
	size_t offset;
	// The slice's length.
	size_t len;
	// The slice's start pointer.
	struct accel_particle *from;
};

// A four-way partition of a 2-dimensional space containing particles.
struct octant;

// An fixed size arena allocator for octants.
struct arena;

// A tree of octants containing particles.
struct particle_tree {
	// The particle tree's root octant.
	struct octant *root;
};

// Recursively constructs the tree structure for the current simulation step.
int particle_tree_build(struct particle_tree *tree,
	const struct accel_particle particles[], float radius, struct arena *arena);

// Executes the current simulation step by updating all particles encompassed
// by the given slice.
//
// Returns the furthest distance to the center of all updated particles.
float particle_tree_simulate(const struct particle_tree *tree,
	const struct particle_slice *slice);

#endif // BARNES_HUT_PHYS_H
