#ifndef __RNG_H__
#define __RNG_H__

#include <stdint.h>
#include <string.h>
#include <stdbool.h>

typedef enum RNG_ERROR
{
	NO_RNG_ERROR,
	ALREADY_INIT_IN_SECURE_MODE,
	RNG_NOT_INITIALIZED,
	UNUSABLE_SEED,
	RNG_FAILED
}
RngError;

typedef enum RNG_MODE
{
	UNKNOWN_MODE,
	NORMAL_MODE,
	SECURE_MODE
}
RngMode;

RngError rngInit(RngMode mode, uint8_t seed[40]);//seed is changing; return new seed

RngError rngGet(uint8_t rng[32]);

//internal using
bool getSystemRandom(uint8_t *rng, size_t length);

#endif /* __RNG_H__ */
