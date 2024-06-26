#include "barnes-hut/options.h"

#include <stdio.h>

#include <getopt.h>

struct options options = {
	.steps	   = 0,
	.particles = 0,
	.max_mass  = 0.0,
	.theta	   = 0.0,
	.threads   = 1,
	.seed	   = 0,
	.optimize  = false,
};

static void print_usage(const char *exe);

int
options_parse(int argc, char *argv[argc])
{
	print_usage(argv[0]);
	return -1;
}

static void
print_usage(const char *exe)
{
	fprintf(stderr,
		"Usage: %s [options]\n"
		"Options:\n"
		"-t [STEPS], --steps=[STEPS]        The number of simulation steps to "
		"compute\n"
		"-n [PARTICLES], --num=[PARTICLES]  The total number of particles "
		"the universe contains.\n"
		"-m [MASS], --mass=[MASS]\n"
		"-r [RADIUS], --radius=[RADIUS]     The initial radius of the "
		"universe.\n"
		"-p [THREADS], --threads=[THREADS]  The number of threads to use for "
		"the simulation\n"
		"-s [SEED], --seed=[SEED]           The seed for random number "
		"generation.\n"
		"-o, --optimize                     xxx"
		"-h, --help                         Print this help and exit\n"
#ifdef DISPLAY
		"-d [DIM], --dimension=[DIM]\n"
#endif // DISPLAY
		"--theta\n",
		exe);
}
