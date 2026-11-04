#ifndef BARNES_HUT_DEF_H
#define BARNES_HUT_DEF_H

#define atomic			_Atomic
#define likely(cond)	__builtin_expect(!!(cond), 1)
#define unlikely(cond)	__builtin_expect(!!(cond), 0)
#define aligned(align)	__attribute__((aligned((align))))
#define printf_like		__attribute__((format(printf, 1, 2)))
#define noreturn		void __attribute__((__noreturn__))
#define unused			__attribute__((unused))

constexpr unsigned cache_line_size = 64;

#define assert(condition, message)
#define abort_oom()		__abort_oom(__FILE__, __func__, __LINE__)

extern noreturn __abort_oom(const char *file, const char *function, int line);

#endif /* BARNES_HUT_DEF_H */
