#ifndef BARNES_HUT_OPTIONS_H
#define BARNES_HUT_OPTIONS_H

#include "barnes-hut/def.h"
#include "barnes-hut/types.h"

struct options {
    /* (t) The number of simulation steps to compute (0 = infinite) */
    unsigned steps;
    /* (n) The total number of particles the universe contains.  */
    unsigned particles;
    /* (m) The mass of each particle */
    double mass;
    /* (r) The initial radius of the universe */
    scalar_t radius;
    /* (p) The number of threads to use for the simulation. */
    unsigned threads;
    /* (s) The seed for random number generation (0..UINT_MAX) */
    unsigned seed;
    /* (o) The step interval for resorting particles (0 = never). */
    unsigned sort_steps;
    /* (f) The flag for enabling a flat universe. */
    bool flat;
    /* (v) The flag for enabling verbose output. */
    bool verbose;
    /* (-) ... */
    scalar_t theta;
    scalar_t dt;
    const char *file;
    bool dump;
    bool no_render;
};

/* The global read-only pointer to the options. */
extern const struct options *options;

extern error_t options_parse(int argc, char *argv[argc]);
extern printf_like void verbose_printf(const char *fmt, ...);

#endif /* BARNES_HUT_OPTIONS_H */
