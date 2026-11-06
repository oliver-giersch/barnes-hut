#include "barnes-hut/options.h"

#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "barnes-hut/types.h"

static struct options options_storage = {
	.steps = 10,
	.verbose = false,
	.flat = false,
};

const struct options *options = &options_storage;

static error_t parse_ull(const char *name, const char *optarg,
	unsigned long long *ull);
static bool parse_ull_failed(const char *optarg, const char *end,
	unsigned long long *ull);
static error_t print_usage(const char binary[]);

extern error_t
options_parse(int argc, char *argv[argc])
{
	enum {
		OPTION_THETA = 1000,
		OPTION_DT,
		OPTION_DUMP,
		OPTION_FILE,
		OPTION_NO_RENDER
	};

	static const char *short_opts = "t:n:m:r:p:s:o:hfv";
	static const char *arg_strs[] = {
		['t'] = "steps",
		['n'] = "num",
		['m'] = "mass",
		['r'] = "radius",
		['p'] = "threads",
		['s'] = "seed",
		['o'] = "optimize",
		[OPTION_THETA] = "theta",
		[OPTION_DT] = "dt",
		[OPTION_NO_RENDER] = "no-render",
	};

	static const struct option long_opts[] = {
		{ "steps", required_argument, NULL, 't' },
		{ "num", required_argument, NULL, 'n' },
		{ "mass", required_argument, NULL, 'm' },
		{ "radius", required_argument, NULL, 'r' },
		{ "threads", required_argument, NULL, 'p' },
		{ "seed", required_argument, NULL, 's' },
		{ "optimize", required_argument, NULL, 'o' },
		{ "help", no_argument, NULL, 'h' },
		{ "flat", no_argument, NULL, 'f' },
		{ "verbose", no_argument, NULL, 'v' },
		{ "theta", required_argument, NULL, OPTION_THETA },
		{ "dt", required_argument, NULL, OPTION_DT },
		{ "file", required_argument, NULL, OPTION_FILE },
		{ "dump", no_argument, NULL, OPTION_DUMP },
		{ "no-render", no_argument, NULL, OPTION_NO_RENDER },
		{ 0, 0, 0, 0 },
	};

	union {
		unsigned long long	ull;
		scalar_t			f;
	} arg;

	error_t error = 0;
	while (true) {
		const int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);

		/* ... */
		if (opt == -1)
			break;

		switch (opt) {
		case 't':
			if ((error = parse_ull(arg_strs[opt], optarg, &arg.ull)))
				goto out;

			if (arg.ull > UINT_MAX) {
				error = EINVAL;
				goto out;
			}

			options_storage.steps = (unsigned)arg.ull;
			break;
		case 'n':
			break;
		case 'm':
			break;
		case 'r':
			break;
		case 'p':
			break;
		}
	}

out:
	if (error)
		print_usage(argv[0]);
	return (error);
}

static error_t
parse_ull(const char *name, const char *optarg, unsigned long long *ull)
{
	if (*optarg == '\0')
		return -1;

	char *end;

	*ull = strtoull(optarg, &end, 10);
	if (!parse_ull_failed(optarg, end, ull))
			return (0);

	error_t error = errno ? errno : EINVAL;
	fprintf(stderr, "invalid argument '%s': '%s' %s\n", name, optarg,
		strerror(error));
	return error;
}

static bool
parse_ull_failed(const char *optarg, const char *end, unsigned long long *ull)
{
	return end == optarg
		|| *end != '\0'
		|| *ull == ULLONG_MAX
		|| errno == ERANGE;
}

static error_t
print_usage(const char binary[])
{
	static constexpr char fmt[] = {
		#embed "usage.text"
		, '\0'
	};

	fprintf(stderr, fmt, binary);
	return ERROR_EXIT_EARLY;
}
