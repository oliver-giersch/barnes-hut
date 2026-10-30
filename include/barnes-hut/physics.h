#ifndef BARNES_HUT_PHYSICS_H
#define BARNES_HUT_PHYSICS_H

#include "barnes-hut/types.h"

/* A 3D-space vector. */
struct vec3 {
	scalar_t x, y, z;
};

/* A point in 3D-space with an associated mass. */
struct point_mass {
	struct vec3 pos;
	mass_t mass;
};

/* A point mass moving in 3D-space with an associated velocity. */
struct particle {
	struct point_mass	pm;
	struct vec3			v;
};

extern size_t particle_count;

extern void particle_init(struct particle *particle, scalar_t x, scalar_t y, scalar_t z, mass_t mass);
extern void particle_randomize(struct particle *particle, scalar_t radius, mass_t mass, bool flat);

#endif /* BARNES_HUT_PHYSICS_H */

