#include "barnes-hut/options.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include "barnes-hut/common.h"

struct options options = {
	.steps	   = 0,
	.particles = 0,
	.max_mass  = 1e12,
	.theta	   = 0.3,
	.threads   = 1,
	.seed	   = 0,
	.optimize  = false,
};

static inline int parse_arg_ull(const char *name, const char *optarg,
	unsigned long long *res);
static inline int parse_arg_float(const char *name, const char *optarg,
	float *res);
static inline bool is_valid_optarg(const char *optarg);

static int print_usage(const char *exe);

// ./barnes $2 $3 256 1e12 0.3 0.01 50 $4 1 1
//         dim,num,r,    m, th, dt, steps,d,o

static const char *argsstrs[] = {
	['t'] = "steps",
	['n'] = "num",
	['m'] = "mass",
	['r'] = "radius",
};

int
options_parse(int argc, char *argv[argc])
{
	static const struct option long_opts[] = {
		{ "steps", required_argument, NULL, 't' },
		{ "num", required_argument, NULL, 'n' },
		{ "mass", required_argument, NULL, 'm' },
		{ "radius", required_argument, NULL, 'r' },
		{ "optimize", no_argument, NULL, 'o' },
		{ 0, 0, 0, 0 },
	};

	int res = 0;
	while (true) {
		const int opt
			= getopt_long(argc, argv, "t:n:m:r:p:s:ho", long_opts, NULL);

		unsigned long long ull;
		// unsigned u;
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
			if (f < 0.0)
				return -1;
			options.max_mass = f;
			break;
		case 'r':
			if ((res = parse_arg_float(argsstrs[opt], optarg, &f)))
				goto out;
			options.radius = f;
			break;
		case 'p':
			ull = strtoull(optarg, NULL, 10);
			break;
		case 's':
			ull = strtoull(optarg, NULL, 10);
			break;
		case 'o':
			options.optimize = true;
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

static inline int
parse_arg_ull(const char *name, const char *optarg, unsigned long long *res)
{
	if (*optarg == '\0')
		return -1;

	char *end;
	*res = strtoull(optarg, &end, 10);
	if (*end != '\0' || *res == ULLONG_MAX) {
		fprintf(stderr, "%s: Invalid %s args: %s\n", __func__, name,
			strerror(errno));
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
	if (end != NULL || *res == HUGE_VALF || *res == -HUGE_VALF) {
		fprintf(stderr, "%s: Invalid %s args: %s\n", __func__, name,
			strerror(errno));
		return errno;
	}

	if (*res < 0.0) {
		fprintf(stderr, "%s: invalid %s args: Negative valute not permitted\n",
			__func__, name);
		return EINVAL;
	}

	return 0;
}

static inline bool
is_valid_optarg(const char *optarg)
{
	return optarg != NULL && *optarg != '\0';
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
		"-o, --optimize                     The flag for enabling memory hierachy optimizations.\n"
		"-h, --help                         Print this help and exit.\n"
#ifdef DISPLAY
		"-d [DIM], --dimension=[DIM]\n"
#endif // DISPLAY
		"--theta\n",
		// clang-format on
		exe);

	return BHE_EARLY_EXIT;
}
