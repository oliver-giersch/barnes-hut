#ifndef BARNES_HUT_OPTIONS_H
#define BARNES_HUT_OPTIONS_H

#include <stdbool.h>
#include <stddef.h>

#include "barnes-hut/common.h"

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
	bool verbose;
} options;

int options_parse(int argc, char *argv[argc]);

printf_like void verbose_printf(const char *fmt, ...);

#endif // BARNES_HUT_OPTIONS_H
