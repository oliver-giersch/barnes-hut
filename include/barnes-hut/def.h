#ifndef BARNES_HUT_DEF_H
#define BARNES_HUT_DEF_H

#define atomic _Atomic
#define likely(cond) __builtin_expect(!!(cond), 1)
#define unlikely(cond) __builtin_expect(!!(cond), 0)
#define aligned(align) __attribute__((aligned((align))))
#define printf_like __attribute__((format(printf, 1, 2)))
#define noreturn void __attribute__((__noreturn__))
#define unused __attribute__((unused))

constexpr unsigned cache_line_size = 64;

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

#endif /* BARNES_HUT_DEF_H */
