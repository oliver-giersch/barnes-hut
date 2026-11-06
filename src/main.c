#include <stdio.h>
#include <stdlib.h>

#include "barnes-hut/def.h"

extern noreturn
__panic(const char *file, const char *function, int line, const char *condition,
	const char *message)
{
	if (message != NULL)
		fprintf(stderr, "%s:%d(%s): assert failed: %s (%s)", function, line,
			file, condition, message);
	else
		fprintf(stderr, "%s:%d(%s): assert failed: %s", function, line, file,
			condition);

	abort();
}