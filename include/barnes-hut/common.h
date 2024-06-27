#ifndef BARNES_HUT_COMMON_H
#define BARNES_HUT_COMMON_H

#define likely(cond) __builtin_expect(!!(cond), 1)
#define unlikely(cond) __builtin_expect(!!(cond), 0)
#define aligned(align) __attribute__((aligned((align))))
#define printf_like __attribute__((format(printf, 1, 2)))

enum error {
	BHE_EARLY_EXIT = -1,
};

#endif // BARNES_HUT_COMMON_H
