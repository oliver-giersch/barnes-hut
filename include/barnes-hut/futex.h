#ifndef BARNES_HUT_FUTEX_H
#define BARNES_HUT_FUTEX_H

#include <stdatomic.h>

static inline error_t futex_wait(atomic_uint32_t *ptr, uint32_t expected);
static inline error_t futex_wake(atomic_uint32_t *ptr, int count);

static inline error_t
futex_wake_all(atomic_uint32_t *ptr)
{
	return futex_wake(ptr, INT_MAX);
}

#endif /* BARNES_HUT_FUTEX_H */

