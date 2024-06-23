#ifndef BH_QUADRANT_H
#define BH_QUADRANT_H

#include <stdbool.h>

struct vec2 {
    float x, y;
};

struct particle {
    struct vec2 pos, vel;
    float mass;
};

struct quadrant {
    struct vec2 center;
    float x, y, len, mass;
    unsigned bodies;
    struct quadrant *children[4];
};

bool quadrant_is_leaf(const struct quadrant *quad);

#endif // BH_QUADRANT_H
