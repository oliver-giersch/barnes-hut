#include "barnes-hut/options.h"

#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include "barnes-hut/common.h"

// The default configuration options.
struct options options = {
	.steps	   = 0,
	.particles = 100000,
	.max_mass  = 1e12,
	.radius	   = 250.0,
	.theta	   = 0.3,
	.dt		   = 0.01,
	.threads   = 1,
	.seed	   = 0,
	.delay	   = 0,
	.optimize  = false,
	.flat	   = false,
	.verbose   = false,
};

static inline int parse_arg_ull(const char *name, const char *optarg,
	unsigned long long *res);
static inline int parse_arg_float(const char *name, const char *optarg,
	float *res);
static int print_usage(const char *exe);

#define THETA 1000
#define DT 1001

static const char *argsstrs[] = {
	['t']	= "steps",
	['n']	= "num",
	['m']	= "mass",
	['r']	= "radius",
	['p']	= "threads",
	['s']	= "seed",
	['d']	= "delay",
	[THETA] = "theta",
	[DT]	= "dt",
};

int
options_parse(int argc, char *argv[argc])
{
	static const struct option long_opts[] = {
		{ "steps", required_argument, NULL, 't' },
		{ "num", required_argument, NULL, 'n' },
		{ "mass", required_argument, NULL, 'm' },
		{ "radius", required_argument, NULL, 'r' },
		{ "theta", required_argument, NULL, THETA },
		{ "dt", required_argument, NULL, DT },
		{ "threads", required_argument, NULL, 'p' },
		{ "seed", required_argument, NULL, 's' },
		{ "delay", required_argument, NULL, 'd' },
		{ "optimize", no_argument, NULL, 'o' },
		{ "flat", no_argument, NULL, 'f' },
		{ "verbose", no_argument, NULL, 'v' },
		{ 0, 0, 0, 0 },
	};

	int res = 0;
	while (true) {
		const int opt
			= getopt_long(argc, argv, "t:n:m:r:p:s:d:hofv", long_opts, NULL);

		unsigned long long ull;
		float f;

		if (opt == -1)
			break;

		switch (opt) {
		case 't':
			if ((res = parse_arg_ull(argsstrs[opt], optarg, &ull)))
				goto out;
			options.steps = (unsigned)ull;
			break;
		case 'n':
			if ((res = parse_arg_ull(argsstrs[opt], optarg, &ull)))
				goto out;
			options.particles = ull;
			break;
		case 'm':
			if ((res = parse_arg_float(argsstrs[opt], optarg, &f)))
				goto out;
			options.max_mass = f;
			break;
		case 'r':
			if ((res = parse_arg_float(argsstrs[opt], optarg, &f)))
				goto out;
			options.radius = f;
			break;
		case 'p':
			if ((res = parse_arg_ull(argsstrs[opt], optarg, &ull)))
				goto out;
			options.threads = (unsigned)ull;
			break;
		case 's':
			if ((res = parse_arg_ull(argsstrs[opt], optarg, &ull)))
				goto out;
			options.seed = (unsigned)ull;
			break;
		case 'd':
			if ((res = parse_arg_ull(argsstrs[opt], optarg, &ull)))
				goto out;
			options.delay = (unsigned)ull;
			break;
		case THETA:
			if ((res = parse_arg_float(argsstrs[opt], optarg, &f)))
				goto out;
			options.theta = f;
			break;
		case DT:
			if ((res = parse_arg_float(argsstrs[opt], optarg, &f)))
				goto out;
			options.dt = f;
			break;
		case 'o':
			options.optimize = true;
			break;
		case 'f':
			options.flat = true;
			break;
		case 'v':
			options.verbose = true;
			break;
		case 'h':
			res = BHE_EARLY_EXIT;
			goto out;
		default:
			break;
		}
	}

out:
	if (res)
		print_usage(argv[0]);
	return res;
}

printf_like void
verbose_printf(const char *fmt, ...)
{
	if (!options.verbose)
		return;

	va_list args;
	va_start(args, fmt);
	(void)vfprintf(stderr, fmt, args);
	va_end(args);
}

static inline int
parse_arg_ull(const char *name, const char *optarg, unsigned long long *res)
{
	if (*optarg == '\0')
		return -1;

	char *end;
	*res = strtoull(optarg, &end, 10);
	if (*end != '\0' || *res == ULLONG_MAX) {
		fprintf(stderr, "Invalid %s args: %s\n", name, strerror(errno));
		return errno;
	}

	return 0;
}

static inline int
parse_arg_float(const char *name, const char *optarg, float *res)
{
	if (*optarg == '\0')
		return -1;

	char *end;
	*res = strtof(optarg, &end);
	if (*end != '\0' || *res == HUGE_VALF || *res == -HUGE_VALF) {
		fprintf(stderr, "Invalid %s arg: %s\n", name, strerror(errno));
		return errno;
	}

	if (*res < 0.0) {
		fprintf(stderr, "Invalid %s arg: Negative valute not permitted\n",
			name);
		return EINVAL;
	}

	return 0;
}

static int
print_usage(const char *exe)
{
	fprintf(stderr,
		// clang-format off
		"Usage: %s [options]\n"
		"Options:\n"
		"-t [STEPS], --steps=[STEPS]        The number of simulation steps to compute.\n"
		"-n [PARTICLES], --num=[PARTICLES]  The total number of particles the universe contains.\n"
		"-m [MASS], --mass=[MASS]           The initial mass of each particle.\n"
		"-r [RADIUS], --radius=[RADIUS]     The initial radius of the universe.\n"
		"-p [THREADS], --threads=[THREADS]  The number of threads to use for the simulation.\n"
		"-s [SEED], --seed=[SEED]           The seed for random number generation (0..UINT_MAX).\n"
		"-d [DELAY], --delay=[DELAY]        The delay in ms after each simulation step.\n"
		"-o, --optimize                     The flag for enabling memory hierachy optimizations.\n"
		"-f, --flat                         The flag for enabling generation of a flat galaxy in the x/y plane.\n"
		"-v, --verbose                      The flag for enabling verbose output.\n"
		"-h, --help                         Print this help and exit.\n"
		"--theta                            The ???\n"
		"--dt                               The g-force dampening factor\n",
		// clang-format on
		exe);

	return BHE_EARLY_EXIT;
}
