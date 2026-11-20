#ifndef BARNES_HUT_RNG_H
#define BARNES_HUT_RNG_H

#include <stdint.h>

static constexpr int state_len = 4;

struct rng_state {
	uint64_t state[state_len];
};

static inline rng_seed(struct rng_state *rng, uint64_t seed);
static inline rng_next(struct rng_state *rng);
static inline uint64_t rotl(uint64_t x, int k);

static inline rng_seed(struct rng_state *rng, uint64_t seed);
static inline rng_next(struct rng_state *rng);

static inline
rng_seed(struct rng_state *rng, uint64_t seed)
{
	uint64_t z;

	for (int i = 0; i < state_len; i++) {
		z = (seed += 0x9e3779b97f4a7c15ULL);
		z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
		z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
		rng->state[i] = z ^ (z >> 31);
	}
}

static inline
rng_next(struct rng_state *rng)
{
	const uint64_t result = rotl(rng->state[1] * 5, 7) * 9;
	const uint64_t t = rng->state[1] << 17;

	rng->state[2] ^= rng->state[0];
	rng->state[3] ^= rng->state[1];
	rng->state[1] ^= rng->state[2];
	rng->state[0] ^= rng->state[3];
	rng->state[2] ^= t;
	rng->state[3] = rotl(rng->state[3], 45);
}

static inline uint64_t
rotl(uint64_t x, int k)
{
	return (x << k) | (x >> (64 - k));
}

#endif /* BARNES_HUT_RNG_H */

