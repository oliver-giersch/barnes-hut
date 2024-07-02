#ifndef BARNES_HUT_ARENA_H
#define BARNES_HUT_ARENA_H

#include <stddef.h>
#include <stdint.h>

#include "barnes-hut/common.h"

#define INVALID_ARENA_ITEM (arena_item_t) UINT32_MAX

typedef uint32_t arena_item_t;

struct arena {
	size_t size;
	size_t item_size;
	arena_item_t curr;
	arena_item_t last;
	void *memory;
};

extern struct arena arena;

static inline void
arena_reset(struct arena *arena)
{
	arena->curr = 0;
}

static inline arena_item_t
arena_malloc(struct arena *arena, size_t size)
{
	arena_item_t item = arena->curr;
	if (unlikely(item == arena->last))
		return INVALID_ARENA_ITEM;

	arena->curr += 1;
	return item;
}

static inline void *
arena_get(struct arena *arena, arena_item_t item)
{
	const uintptr_t addr
		= ((uintptr_t)arena->memory + (uintptr_t)item * arena->item_size);
	return (void *)addr;
}

int arena_init(struct arena *arena, size_t size, size_t item_size);
void arena_deinit(struct arena *arena);

#endif // BARNES_HUT_ARENA_H
