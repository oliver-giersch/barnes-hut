#ifndef BARNES_HUT_ARENA_H
#define BARNES_HUT_ARENA_H

#include <stddef.h>
#include <stdint.h>

#include "../include/common.h"

struct arena {
	void *memory;
	void *curr;
	void *end;
};

static inline void
arena_reset(struct arena *arena)
{
	arena->curr = arena->memory;
}

static inline void *
arena_malloc(struct arena *arena, size_t size)
{
	void *item = arena->curr;
	if (unlikely(item > arena->end))
		return NULL;
	arena->curr = (void *)((uintptr_t)arena->curr + size);
	return item;
}

int arena_init(struct arena *arena, size_t size);
void arena_deinit(struct arena *arena);

#endif // BARNES_HUT_ARENA_H
