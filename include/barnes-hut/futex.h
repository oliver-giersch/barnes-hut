#ifndef BARNES_HUT_FUTEX_H
#define BARNES_HUT_FUTEX_H

#include <limits.h>
#include <stdatomic.h>
#include <stdint.h>

#include "barnes-hut/def.h"
#include "barnes-hut/types.h"

typedef atomic uint32_t atomic_uint32_t;

static inline error_t futex_wait(atomic_uint32_t *ptr, uint32_t expected);
static inline error_t futex_wake(atomic_uint32_t *ptr, int count);

static inline error_t
futex_wake_all(atomic_uint32_t *ptr)
{
	return futex_wake(ptr, INT_MAX);
}

#if defined(__linux__)
#include "barnes-hut/futex/linux.h"
#elif defined(__FreeBSD__)
#include "barnes-hut/futex/freebsd.h"
#else
#error "unsupported OS"
#endif /* __linux__ || __FreeBSD__ */

#endif /* BARNES_HUT_FUTEX_H */

