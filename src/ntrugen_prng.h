#ifndef NTRUGEN_PRNG_H__
#define NTRUGEN_PRNG_H__

#include <stddef.h>
#include <stdint.h>

#include "sha3.h"

/*
 * This file implements a PRNG based on ChaCha8, as well as an API to
 * access the system RNG.
 */

/*
 * Get a high-entropy random seed from the operating system. This function
 * return 0 on success (all requested bytes have been obtained), or a
 * negative error code on error.
 */
int ntrugen_sysrng(void *seed, size_t len);

/*
 * State structure for a ChaCha8-based RNG. Contents are opaque.
 */
typedef struct {
	union {
		uint8_t d[40];
		uint64_t dummy_u64;
	} state;
} ntrugen_prng_chacha8_context;

/*
 * Initialize a ChaCha8 PRNG. The provided seed should contain enough
 * entropy for the expected security target. The seed length MUST NOT
 * exceed 32 bytes.
 */
void ntrugen_prng_chacha8_init(ntrugen_prng_chacha8_context *ctx,
	const void *seed, size_t seed_len);

/*
 * Generate some more pseudorandom bytes.
 *
 * The 'ctx' parameter MUST point to a ntrugen_prng_chacha8_context structure;
 * it is here declared under type 'void *' so as to be compatible with
 * the generic ntrugen_rng type (declared in ntrugen.h).
 *
 * This implementation internally generates bytes by relatively large chunks
 * (64 or 512 bytes, depending on the architecture); calling this function
 * repeatedly for small outputs is inefficient.
 */
void ntrugen_prng_chacha8_out(void *ctx, void *dst, size_t len);

/*
 * State structure for a SHAKE-based PRNG. Four SHAKE256 contexts are used
 * internally, with interleaved outputs (on a 64-bit granularity).
 */
typedef struct {
	shake_x4_context state;
} ntrugen_prng_shake_context;

/*
 * Initialize a SHAKE256 PRNG. The provided seed should contain enough
 * entropy for the expected security target. There is no practical upper
 * limit on seed length, but seed sizes of at most 134 bytes are best
 * for performance.
 */
void ntrugen_prng_shake_init(ntrugen_prng_shake_context *ctx,
	const void *seed, size_t seed_len);

/*
 * Generate some more pseudorandom bytes.
 *
 * The 'ctx' parameter MUST point to a ntrugen_prng_shake_context structure;
 * it is here declared under type 'void *' so as to be compatible with
 * the generic ntrugen_rng type (declared in ntrugen.h).
 *
 * This implementation internally generates bytes by chunks of 32 bytes;
 * calling this function repeatedly for small outputs is inefficient.
 */
void ntrugen_prng_shake_out(void *ctx, void *dst, size_t len);

#endif
