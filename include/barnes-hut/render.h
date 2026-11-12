#ifndef BARNES_HUT_RENDER_H
#define BARNES_HUT_RENDER_H

#include "barnes-hut/types.h"

#ifdef OPTION_RENDER
extern error_t render_main(void); // fixme: arguments?

extern error_t render_init(void);
extern void render_destroy(void);
extern bool render_scene(const struct particle particles[], scalar_t radius);
extern bool render_wait(void);
extern bool render_wake(void);
#else
static inline error_t
render_init(void)
{
	return 0;
}

static inline void
render_destroy(void) {}

static inline bool
render_scene(const struct particle particles[], scalar_t radius)
{
	return true;
}
#endif /* OPTION_RENDER */

#endif /* BARNES_HUT_RENDER_H */
