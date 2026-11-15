#ifndef BARNES_HUT_RENDER_H
#define BARNES_HUT_RENDER_H

#include "barnes-hut/types.h"

#ifdef OPTION_RENDER
extern error_t render_main(void); // fixme: arguments?

extern bool render_wait(void);
extern bool render_wake(scalar_t radius);
#else
static inline bool
render_wait(void)
{
	return true;
}

static inline bool
render_wake(scalar_t radius)
{
	(void)radius;
	return true;
}
#endif /* OPTION_RENDER */
#endif /* BARNES_HUT_RENDER_H */
