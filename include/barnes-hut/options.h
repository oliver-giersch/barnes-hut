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
	float theta;
	unsigned threads;
	unsigned seed;
	bool optimize;
	bool verbose;
} options;

int options_parse(int argc, char *argv[argc]);

printf_like void verbose_printf(const char *fmt, ...);

#endif // BARNES_HUT_OPTIONS_H
