#ifndef BARNES_HUT_TYPES_H
#define BARNES_HUT_TYPES_H

/* The type representing a coordinate in space. */
#if (OPTION_DOUBLE_PRECISION == 1)
typedef float scalar_t;
#else
typedef double scalar_t;
#endif /* OPTION_DOUBLE_PRECISION */

/* The type describing a particle's mass. */
typedef double mass_t;

struct vec3;
struct particle;
struct octant;

#endif /* BARNES_HUT_TYPES_H */

