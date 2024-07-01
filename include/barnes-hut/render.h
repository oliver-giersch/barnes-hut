#ifndef BARNES_HUT_RENDER_H
#define BARNES_HUT_RENDER_H

#include <stdbool.h>

#include "barnes-hut/phys.h"

int render_init(void);
void render_deinit(void);
bool render_scene(const struct accel_particle particles[], float radius);

#endif // BARNES_HUT_RENDER_H
