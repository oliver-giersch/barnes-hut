#ifndef BARNES_HUT_FUTEX_LINUX_H
#define BARNES_HUT_FUTEX_LINUX_H

#include <sys/syscall.h>
#include <linux/futex.h>

#include "barnes-hut/def.h"

static inline error_t
futex_wait(atomic_uint32_t *ptr, uint32_t expected)
{
	const long res = syscall(SYS_FUTEX, ptr, FUTEX_WAIT, expected, nullptr);
	return (error_t)res;
}

static inline error_t
futex_wake(atomic_uint32_t *ptr, int count)
{
	const long res = syscall(SYS_FUTEX, ptr, FUTEX_WAKE, (uint32_t)count,
		nullptr);
	return (error_t)res;
}

#endif /* BARNES_HUT_FUTEX_LINUX_H */
