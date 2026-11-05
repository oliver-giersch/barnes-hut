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

/* The standard error type. */
typedef int error_t;

enum : error_t {
    ERROR_SUCCESS		= 0,
    ERROR_EXIT_EARLY	= 1,
};

struct vec3;
struct particle;
struct octant;
struct octant_cache;

#endif /* BARNES_HUT_TYPES_H */
