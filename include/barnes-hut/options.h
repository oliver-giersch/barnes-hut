#ifndef BARNES_HUT_OPTIONS_H
#define BARNES_HUT_OPTIONS_H

#include "barnes-hut/def.h"
#include "barnes-hut/types.h"

struct options {
    /* The number of discrete steps to simulate (0 = infinity). */
    unsigned steps;
    /* The flag for enabling verbose output. */
    bool verbose;
    /* The flag for enabling a flat universe. */
    bool flat;
};

/* The global read-only pointer to the options. */
extern const struct options *options;

extern error_t options_parse(int argc, char *argv[argc]);
extern printf_like void verbose_printf(const char *fmt, ...);

#endif /* BARNES_HUT_OPTIONS_H */
