#ifndef MT19937_64_H
#define MT19937_64_H

// Initializes RNG with a seed.
void mt1993764_init(unsigned long long seed);
// Generates a random number in an [0, 2^64-1] interval.
unsigned long long mt1993764_int64(void);

#endif // MT19937_64_H
