/*
 * PRNG and interface to the system RNG.
 */

#include "ntrugen_prng.h"
#include "ng_inner.h"

/*
 * Include relevant system header files. For Win32, this will also need
 * linking with advapi32.dll, which we trigger with an appropriate #pragma.
 */
#if NTRUGEN_RAND_GETENTROPY
#include <unistd.h>
#endif
#if NTRUGEN_RAND_URANDOM
#include <sys/types.h>
#if !NTRUGEN_RAND_GETENTROPY
#include <unistd.h>
#endif
#include <fcntl.h>
#include <errno.h>
#endif
#if NTRUGEN_RAND_WIN32
#include <windows.h>
#define SystemFunction036   NTAPI SystemFunction036
#include <NTSecAPI.h>
#undef SystemFunction036
#pragma comment(lib, "advapi32")
#endif

/* see ntrugen_prng.h */
int
ntrugen_sysrng(void *seed, size_t len)
{
	(void)seed;
	if (len == 0) {
		return 0;
	}
#if NTRUGEN_RAND_GETENTROPY
	while (len > 0) {
		size_t clen = len < 256 ? len : 256;
		if (getentropy(seed, clen) != 0) {
			break;
		}
		seed = (uint8_t *)seed + clen;
		len -= clen;
	}
	if (len == 0) {
		return 0;
	}
#endif
#if NTRUGEN_RAND_URANDOM
	int f = open("/dev/urandom", O_RDONLY);
	if (f >= 0) {
		while (len > 0) {
			ssize_t rlen = read(f, seed, len);
			if (rlen < 0) {
				if (errno == EINTR) {
					continue;
				}
				break;
			}
			seed = (uint8_t *)seed + rlen;
			len -= (size_t)rlen;
		}
		close(f);
		if (len == 0) {
			return 0;
		}
	}
#endif
#if NTRUGEN_RAND_WIN32
	while (len > 0) {
		size_t clen = len < 256 ? len : 256;
		if (!RtlGenRandom(seed, clen)) {
			break;
		}
		seed = (uint8_t *)seed + clen;
		len -= clen;
	}
	if (len == 0) {
		return 0;
	}
	/*
	 * The nominal Win32 method is to use CryptoAPI, but RtlGenRandom()
	 * is more efficient and works since Windows XP and Server 2003.
	HCRYPTPROV hp;
	if (CryptAcquireContext(&hp, 0, 0, PROV_RSA_FULL,
		CRYPT_VERIFYCONTEXT | CRYPT_SILENT))
	{
		BOOL r = CryptGenRandom(hp, (DWORD)len, seed);
		CryptReleaseContext(hp, 0);
		if (r) {
			return 0;
		}
	}
	 */
#endif
	return -1;
}

/* see ntrugen_prng.h */
void
ntrugen_prng_chacha8_init(ntrugen_prng_chacha8_context *ctx,
	const void *seed, size_t seed_len)
{
	/*
	 * Truncate the seed to 32 bytes; it is incorrect to provide a larger
	 * seed, but we'd prefer not having a buffer overflow in that case.
	 */
	if (seed_len > 32) {
		seed_len = 32;
	}

	/*
	 * We use the seed as a key in ChaCha8 (padded with zeros if it is
	 * shorter than 32 bytes). To distinguish it from ChaCha20 as used
	 * for encryption purposes (RFC 8439) (beyond the reduced number
	 * of rounds), we modify the first constant word.
	 *
	 * To ensure reproducibility for a given seed, we enforce
	 * little-endian interpretation of the seed words.
	 */
	uint8_t tmp[32];
	memcpy(tmp, seed, seed_len);
	memset(tmp + seed_len, 0, (sizeof tmp) - seed_len);
	uint32_t *key = (uint32_t *)&ctx->state;
	for (size_t u = 0; u < 8; u ++) {
		key[u] = dec32le(tmp + (u << 2));
	}
	*(uint64_t *)(key + 8) = 0;
}

/* see ntrugen_prng.h */
TARGET_AVX2
void
ntrugen_prng_chacha8_out(void *ctx, void *dst, size_t len)
{
	/*
	 * The AVX2 and portable implementations agree with each other
	 * as long as all calls are for len == 512. If a call is made
	 * with a different length then the two implementations diverge
	 * (though both are safe as a PRNG). The ntrugen code only invokes
	 * the PRNG with len == 512.
	 */

	static const uint32_t CW[] = {
		/* Note: first constant differs from RFC 8439, in case
		   somebody is misguided enough to use as seed the same
		   value as a key used for encryption. */
		0xA7C083FE, 0x3320646E, 0x79622d32, 0x6B206574
	};

	uint8_t *buf = dst;
	uint64_t cc = *((uint64_t *)ctx + 4);

#if NTRUGEN_AVX2
	while (len > 0) {
		__m256i state[16], init[16];
		uint32_t *sw = (uint32_t *)ctx;

		for (size_t u = 0; u < 4; u ++) {
			state[u] = init[u] = _mm256_set1_epi32(CW[u]);
		}
		for (size_t u = 4; u < 12; u ++) {
			state[u] = init[u] = _mm256_set1_epi32(sw[u - 4]);
		}
		state[12] = init[12] = _mm256_setr_epi32(
			(uint32_t)(cc + 0), (uint32_t)(cc + 1),
			(uint32_t)(cc + 2), (uint32_t)(cc + 3),
			(uint32_t)(cc + 4), (uint32_t)(cc + 5),
			(uint32_t)(cc + 6), (uint32_t)(cc + 7));
		state[13] = init[13] = _mm256_setr_epi32(
			(uint32_t)((cc + 0) >> 32), (uint32_t)((cc + 1) >> 32),
			(uint32_t)((cc + 2) >> 32), (uint32_t)((cc + 3) >> 32),
			(uint32_t)((cc + 4) >> 32), (uint32_t)((cc + 5) >> 32),
			(uint32_t)((cc + 6) >> 32), (uint32_t)((cc + 7) >> 32));
		cc += 8;
		state[14] = init[14] = _mm256_setzero_si256();
		state[15] = init[15] = _mm256_setzero_si256();

		/*
		 * Do all rounds.
		 */
		for (int i = 0; i < 4; i ++) {

#define QROUND(a, b, c, d)   do { \
		state[a] = _mm256_add_epi32(state[a], state[b]); \
		state[d] = _mm256_xor_si256(state[d], state[a]); \
		state[d] = _mm256_or_si256( \
			_mm256_slli_epi32(state[d], 16), \
			_mm256_srli_epi32(state[d], 16)); \
		state[c] = _mm256_add_epi32(state[c], state[d]); \
		state[b] = _mm256_xor_si256(state[b], state[c]); \
		state[b] = _mm256_or_si256( \
			_mm256_slli_epi32(state[b], 12), \
			_mm256_srli_epi32(state[b], 20)); \
		state[a] = _mm256_add_epi32(state[a], state[b]); \
		state[d] = _mm256_xor_si256(state[d], state[a]); \
		state[d] = _mm256_or_si256( \
			_mm256_slli_epi32(state[d],  8), \
			_mm256_srli_epi32(state[d], 24)); \
		state[c] = _mm256_add_epi32(state[c], state[d]); \
		state[b] = _mm256_xor_si256(state[b], state[c]); \
		state[b] = _mm256_or_si256( \
			_mm256_slli_epi32(state[b], 7), \
			_mm256_srli_epi32(state[b], 25)); \
	} while (0)

			QROUND( 0,  4,  8, 12);
			QROUND( 1,  5,  9, 13);
			QROUND( 2,  6, 10, 14);
			QROUND( 3,  7, 11, 15);
			QROUND( 0,  5, 10, 15);
			QROUND( 1,  6, 11, 12);
			QROUND( 2,  7,  8, 13);
			QROUND( 3,  4,  9, 14);

#undef QROUND

		}

		/*
		 * Add initial state back and encode the result in the
		 * destination buffer. We can dump the AVX2 values "as
		 * is" because the non-AVX2 code uses a compatible order
		 * of values.
		 */
		if (len >= 512) {
			for (size_t u = 0; u < 16; u ++) {
				_mm256_storeu_si256((__m256i *)&buf[u << 5],
					_mm256_add_epi32(state[u], init[u]));
			}
			buf += 512;
			len -= 512;
		} else {
			for (size_t u = 0; len > 0; u ++) {
				union {
					uint8_t b[32];
					__m256i y;
				} t;
				t.y = _mm256_add_epi32(state[u], init[u]);
				if (len >= 32) {
					memcpy(buf, t.b, 32);
					buf += 32;
					len -= 32;
				} else {
					memcpy(buf, t.b, len);
					len = 0;
				}
			}
			break;
		}
	}

#else
	/*
	 * Portable version.
	 */

	/*
	 * This is normally invoked by chunks of 512 bytes; in that case,
	 * we interleave the words of the eight successive 64-byte outputs
	 * so as to get the same final output as the AVX2 implementation.
	 */
	int out512 = (len == 512);

	/*
	 * State uses local endianness. Only the output bytes must be
	 * converted to little endian (if used on a big-endian machine).
	 */
	while (len > 0) {
		uint32_t state[16];
		memcpy(&state[0], CW, sizeof CW);
		memcpy(&state[4], ctx, 32);
		state[12] = (uint32_t)cc;
		state[13] = (uint32_t)(cc >> 32);
		state[14] = 0;
		state[15] = 0;
		for (int i = 0; i < 4; i ++) {
#define QROUND(a, b, c, d)   do { \
		state[a] += state[b]; \
		state[d] ^= state[a]; \
		state[d] = (state[d] << 16) | (state[d] >> 16); \
		state[c] += state[d]; \
		state[b] ^= state[c]; \
		state[b] = (state[b] << 12) | (state[b] >> 20); \
		state[a] += state[b]; \
		state[d] ^= state[a]; \
		state[d] = (state[d] <<  8) | (state[d] >> 24); \
		state[c] += state[d]; \
		state[b] ^= state[c]; \
		state[b] = (state[b] <<  7) | (state[b] >> 25); \
	} while (0)

			QROUND( 0,  4,  8, 12);
			QROUND( 1,  5,  9, 13);
			QROUND( 2,  6, 10, 14);
			QROUND( 3,  7, 11, 15);
			QROUND( 0,  5, 10, 15);
			QROUND( 1,  6, 11, 12);
			QROUND( 2,  7,  8, 13);
			QROUND( 3,  4,  9, 14);

#undef QROUND
		}

		for (size_t v = 0; v < 4; v ++) {
			state[v] += CW[v];
		}
		for (size_t v = 4; v < 12; v ++) {
			state[v] += ((uint32_t *)ctx)[v - 4];
		}
		state[12] += (uint32_t)cc;
		state[13] += (uint32_t)(cc >> 32);
		cc ++;

		if (out512) {
			for (size_t v = 0; v < 16; v ++) {
				enc32le(buf + (v << 5), state[v]);
			}
			buf += 4;
			len -= 64;
		} else {
			if (len >= 64) {
				for (size_t v = 0; v < 16; v ++) {
					enc32le(buf + (v << 2), state[v]);
				}
				buf += 64;
				len -= 64;
			} else {
				size_t v;
				for (v = 0; len >= 4;
					v ++, buf += 4, len -= 4)
				{
					enc32le(buf, state[v]);
				}
				uint32_t x = state[v];
				while (len > 0) {
					*buf ++ = (uint8_t)x;
					x >>= 8;
					len --;
				}
				break;
			}
		}
	}

#endif

	*((uint64_t *)ctx + 4) = cc;
}

/* see ntrugen_prng.h */
void
ntrugen_prng_shake_init(ntrugen_prng_shake_context *ctx,
	const void *seed, size_t seed_len)
{
	shake_context sc[4];
	shake_init(&sc[0], 256);
	shake_inject(&sc[0], seed, seed_len);
	for (int i = 1; i < 4; i ++) {
		sc[i] = sc[0];
	}
	for (int i = 0; i < 4; i ++) {
		uint8_t x = i;
		shake_inject(&sc[i], &x, 1);
	}
	shake_x4_flip(&ctx->state, sc);
}

/* see ntrugen_prng.h */
void
ntrugen_prng_shake_out(void *ctx, void *dst, size_t len)
{
	shake_x4_context *scx4 = &((ntrugen_prng_shake_context *)ctx)->state;
	uint8_t *buf = dst;
#if NTRUGEN_LE && NTRUGEN_UNALIGNED
	shake_x4_extract_words(scx4, dst, len >> 5);
	buf += (len & ~(size_t)31);
	len &= 31;
#endif
	while (len > 0) {
		uint64_t tq[4];
		shake_x4_extract_words(scx4, tq, 1);
		size_t clen = len < 32 ? len : 32;
		for (size_t u = 0; u < clen; u ++) {
			buf[u] = (uint8_t)(tq[u >> 3] >> ((u & 7) << 3));
		}
		buf += clen;
		len -= clen;
	}
}
