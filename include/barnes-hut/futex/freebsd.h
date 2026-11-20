#ifndef BARNES_HUT_FUTEX_FREEBSD_H
#define BARNES_HUT_FUTEX_FREEBSD_H

#include <sys/types.>
#include <sys/umtx.h>

#include "barnes-hut/def.h"

static inline error_t
futex_wait(atomic_uint32_t *ptr, uint32_t expected)
{
	const int res = _umtx_op((void *)ptr, UMTX_OP_WAIT_UINT_PRIVATE, expected,
		nullptr, nullptr);
	// handle some errors?
	return res;
}

static inline error_t
futex_wake(atomic_uint32_t *ptr, int count)
{
	const int res = _umtx_op((void *)ptr, UMTX_OP_WAKE_PRIVATE, count, nullptr,
		nullptr);
	// handle some errors?
	return res;
}

#endif /* BARNES_HUT_FUTEX_FREEBSD_H */

