#ifndef BARNES_HUT_DEF_H
#define BARNES_HUT_DEF_H

#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif /* __x86_64__ || __i386__ */

/* Attribute and built-in wrappers. */

#define atomic _Atomic
#define likely(cond) __builtin_expect(!!(cond), 1)
#define unlikely(cond) __builtin_expect(!!(cond), 0)
#define aligned(align) __attribute__((aligned((align))))
#define printf_like __attribute__((format(printf, 1, 2)))
#define noreturn void __attribute__((__noreturn__))
#define unused __attribute__((unused))

/* Build options */

#ifdef OPTION_DOUBLE_PRECISION
	static constexpr bool double_precision = true;
#else
	static constexpr bool double_precision = false;
#endif /* OPTION_DOUBLE_PRECISION */

#ifdef OPTION_LARGE_SIMULATION
	static constexpr bool large_simulation = true;
#else
	static constexpr bool large_simulation = false;
#endif /* OPTION_LARGE_SIMULATION */

#ifdef OPTION_NUMA_AWARE
	static constexpr bool numa_aware = true;
#else
	static constexpr bool numa_aware = false;
#endif /* OPTION_NUMA_AWARE */

#ifdef OPTION_PRINT_PROGRESS
	static constexpr bool print_progress = true;
#else
	static constexpr bool print_progress = false;
#endif /* OPTION_PRINT_PROGRESS */

#ifdef OPTION_RENDER
	static constexpr bool render = true;
#else
	static constexpr bool render = false;
#endif /* OPTION_RENDER */

static constexpr unsigned cache_line_size = 64;

/* Aborts the process on an OOM condition. */
#define abort_oom() __abort_oom(__FILE__, __func__, __LINE__)

/*
 * Asserts that condition is true and aborts the process otherwise after
 * printing the given message to stderr.
 */
#ifndef NDEBUG
#define assert(cond, message)                                                  \
	do {                                                                       \
		if (unlikely(!(cond)))                                                 \
			__panic(__FILE__, __func__, __LINE__, (#cond), (message));         \
	} while (false)
#else
#define assert(cond, message) ((void)0)
#endif /* NDEBUG */

extern noreturn __panic(const char *file, const char *function, int line,
	const char *condition, const char *message);
extern noreturn __abort_oom(const char *file, const char *function, int line);

static inline void
cpu_relax(void)
{
#if defined(__x86_64__) || defined(__i386__)
	_mm_pause();
#endif /* __x86_64__ || __i386__ */
}

#endif /* BARNES_HUT_DEF_H */
