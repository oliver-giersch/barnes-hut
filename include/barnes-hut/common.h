#ifndef BARNES_HUT_COMMON_H
#define BARNES_HUT_COMMON_H

#define likely(cond) __builtin_expect(!!(cond), 1)
#define unlikely(cond) __builtin_expect(!!(cond), 0)
#define aligned(align) __attribute__((aligned((align))))

enum error {
	ENOMEM = -1,
};

#endif // BARNES_HUT_COMMON_H
