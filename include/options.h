#ifndef BARNES_HUT_OPTIONS_H
#define BARNES_HUT_OPTIONS_H

#include <stddef.h>

// The global options and settings.
extern struct options {
	unsigned steps;
	size_t particles;
	float max_mass;
	float radius;
	float theta;
	unsigned threads;
	unsigned seed;
	bool optimize;
} options;

int options_parse(int argc, char *argv[argc]);

#endif // BARNES_HUT_OPTIONS_H
