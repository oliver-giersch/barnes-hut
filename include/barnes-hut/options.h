#ifndef BARNES_HUT_OPTIONS_H
#define BARNES_HUT_OPTIONS_H

#include <stdbool.h>
#include <stddef.h>

#include "barnes-hut/common.h"

// The global options and settings.
extern struct options {
	// The number of simulation steps to perform (0 means infinite).
	unsigned steps;
	// The number of particles to simulate.
	size_t particles;
	// The initial mass of each particle.
	float max_mass;
	// The initial radius of the galaxy.
	float radius;
	// The ???.
	float theta;
	// The total number of threads to utilize.
	unsigned threads;
	// The seed for RNG (0 means no fixed seed).
	unsigned seed;
	// The delay in ms afer each simulation step.
	unsigned delay;
	// The flag for enabling z-curve order sorting optimization.
	bool optimize;
	// The flag for enabling more verbose output to `stderr`.
	bool verbose;
	// The flag for forcing the entire galaxy into a flat x/y plane.
	bool flat;
} options;

int options_parse(int argc, char *argv[argc]);

printf_like void verbose_printf(const char *fmt, ...);

#endif // BARNES_HUT_OPTIONS_H
