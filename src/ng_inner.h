#ifndef NG_INNER_H__
#define NG_INNER_H__

/* ==================================================================== */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "ng_config.h"

#ifndef NTRUGEN_PREFIX
#define NTRUGEN_PREFIX   ntrugen
#endif
#define Zn(name)             Zn_(NTRUGEN_PREFIX, name)
#define Zn_(prefix, name)    Zn__(prefix, name)
#define Zn__(prefix, name)   prefix ## _ ## name

#include "sha3.h"

/* ==================================================================== */

#ifndef NTRUGEN_AVX2
/*
 * Auto-detection of AVX2 support, if not overridden by configuration.
 */
#if defined __AVX2__ && __AVX2__
#define NTRUGEN_AVX2   1
#else
#define NTRUGEN_AVX2   0
#endif
#endif // NTRUGEN_AVX2

#if NTRUGEN_AVX2
/*
 * This implementation uses AVX2 intrinsics.
 */
#include <immintrin.h>
#if defined __GNUC__ || defined __clang__
#include <x86intrin.h>
#endif
#ifndef NTRUGEN_LE
#define NTRUGEN_LE   1
#endif
#ifndef NTRUGEN_UNALIGNED
#define NTRUGEN_UNALIGNED   1
#endif
#if defined __GNUC__
#define TARGET_AVX2    __attribute__((target("avx2,lzcnt,pclmul")))
#define ALIGNED_AVX2   __attribute__((aligned(32)))
#elif defined _MSC_VER && _MSC_VER
#pragma warning( disable : 4752 )
#endif
#endif // NTRUGEN_AVX2

#ifndef TARGET_AVX2
#define TARGET_AVX2
#endif
#ifndef ALIGNED_AVX2
#define ALIGNED_AVX2
#endif

/*
 * Auto-detect ARM Cortex-M4 platforms.
 */
#ifndef NTRUGEN_ASM_CORTEXM4
#if (defined __ARM_ARCH_7EM__ && __ARM_ARCH_7EM__) \
	&& (defined __ARM_FEATURE_DSP && __ARM_FEATURE_DSP)
#define NTRUGEN_ASM_CORTEXM4   1
#else
#define NTRUGEN_ASM_CORTEXM4   0
#endif
#endif

/*
 * Disable warning on applying unary minus on an unsigned type.
 */
#if defined _MSC_VER && _MSC_VER
#pragma warning( disable : 4146 )
#pragma warning( disable : 4244 )
#pragma warning( disable : 4267 )
#pragma warning( disable : 4334 )
#endif

/*
 * Auto-detect 64-bit architectures.
 */
#ifndef NTRUGEN_64
#if defined __x86_64__ || defined _M_X64 \
        || defined __ia64 || defined __itanium__ || defined _M_IA64 \
        || defined __powerpc64__ || defined __ppc64__ || defined __PPC64__ \
        || defined __64BIT__ || defined _LP64 || defined __LP64__ \
        || defined __sparc64__ \
        || defined __aarch64__ || defined _M_ARM64 \
        || defined __mips64
#define NTRUGEN_64   1
#else
#define NTRUGEN_64   0
#endif
#endif

/*
 * Auto-detect endianness and support of unaligned accesses.
 */
#if defined __i386__ || defined _M_IX86 \
        || defined __x86_64__ || defined _M_X64 \
	|| defined __aarch64__ || defined _M_ARM64 || defined _M_ARM64EC \
        || (defined _ARCH_PWR8 \
                && (defined __LITTLE_ENDIAN || defined __LITTLE_ENDIAN__))

#ifndef NTRUGEN_LE
#define NTRUGEN_LE   1
#endif
#ifndef NTRUGEN_UNALIGNED
#define NTRUGEN_UNALIGNED   1
#endif

#elif (defined __LITTLE_ENDIAN || defined __LITTLE_ENDIAN__) \
        || (defined __BYTE_ORDER__ && defined __ORDER_LITTLE_ENDIAN__ \
                && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)

#ifndef NTRUGEN_LE
#define NTRUGEN_LE   1
#endif
#ifndef NTRUGEN_UNALIGNED
#define NTRUGEN_UNALIGNED   0
#endif

#else

#ifndef NTRUGEN_LE
#define NTRUGEN_LE   0
#endif
#ifndef NTRUGEN_UNALIGNED
#define NTRUGEN_UNALIGNED   0
#endif

#endif

/*
 * For seed generation:
 *
 *  - On Linux (glibc-2.25+), FreeBSD 12+ and OpenBSD, use getentropy().
 *  - On other Unix-like systems, use /dev/urandom (also a fallback for
 *    failed getentropy() calls).
 *  - On Windows, use CryptGenRandom().
 */

#ifndef NTRUGEN_RAND_GETENTROPY
#if (defined __linux && defined __GLIBC__ \
        && (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 25))) \
        || (defined __FreeBSD__ && __FreeBSD__ >= 12) \
        || defined __OpenBSD__
#define NTRUGEN_RAND_GETENTROPY   1
#else
#define NTRUGEN_RAND_GETENTROPY   0
#endif
#endif
#ifndef NTRUGEN_RAND_URANDOM
#if defined _AIX \
        || defined __ANDROID__ \
        || defined __FreeBSD__ \
        || defined __NetBSD__ \
        || defined __OpenBSD__ \
        || defined __DragonFly__ \
        || defined __linux__ \
        || (defined __sun && (defined __SVR4 || defined __svr4__)) \
        || (defined __APPLE__ && defined __MACH__)
#define NTRUGEN_RAND_URANDOM   1
#else
#define NTRUGEN_RAND_URANDOM   0
#endif
#endif

#ifndef NTRUGEN_RAND_WIN32
#if defined _WIN32 || defined _WIN64
#define NTRUGEN_RAND_WIN32   1
#else
#define NTRUGEN_RAND_WIN32   0
#endif
#endif

/*
 * MSVC 2015 does not known the C99 keyword 'restrict'.
 */
#if defined _MSC_VER && _MSC_VER
#ifndef restrict
#define restrict   __restrict
#endif
#endif

/*
 * Enable stats if not defined explicitly (stats should be disabled by
 * default, since they are not thread-safe).
 */
#ifndef NTRUGEN_STATS
#define NTRUGEN_STATS   0
#endif

/* ==================================================================== */
/*
 * The ntrugen public API is redeclared here, but with the proper name
 * prefixes and macros.
 */

typedef void (*ntrugen_rng)(void *ctx, void *dst, size_t len);

#define BAT_keygen   Zn(BAT_keygen)
int BAT_keygen(unsigned logn,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G, int32_t *w,
	ntrugen_rng rng, void *rng_context,
	void *tmp, size_t tmp_len);

#define BAT_recover_G   Zn(BAT_recover_G)
int BAT_recover_G(unsigned logn,
	int8_t *G,
	const int8_t *f, const int8_t *g, const int8_t *F,
	void *tmp, size_t tmp_len);

#define Falcon_keygen   Zn(Falcon_keygen)
int Falcon_keygen(unsigned logn,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G,
	ntrugen_rng rng, void *rng_context,
	void *tmp, size_t tmp_len);

#define Hawk_keygen   Zn(Hawk_keygen)
int Hawk_keygen(unsigned logn,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G,
	int16_t *q00, int16_t *q01, int32_t *q11,
	void *seed, ntrugen_rng rng, void *rng_context,
	void *tmp, size_t tmp_len);

#define Hawk_regen_fg   Zn(Hawk_regen_fg)
void Hawk_regen_fg(unsigned logn,
	int8_t *f, int8_t *g, const void *seed);

#define Hawk_recover_G   Zn(Hawk_recover_G)
int Hawk_recover_G(unsigned logn,
	int8_t *G,
	const int8_t *f, const int8_t *g, const int8_t *F,
	void *tmp, size_t tmp_len);

#define Hawk_recover_qq   Zn(Hawk_recover_qq)
int Hawk_recover_qq(unsigned logn,
	int16_t *q00, int16_t *q01, int32_t *q11,
	const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G,
	void *tmp, size_t tmp_len);

/* ==================================================================== */

static inline unsigned
dec16le(const void *src)
{
#if NTRUGEN_LE && NTRUGEN_UNALIGNED
	return *(const uint16_t *)src;
#else
	const uint8_t *buf = src;
	return (unsigned)buf[0]
		| ((unsigned)buf[1] << 8);
#endif
}

static inline void
enc16le(void *dst, unsigned x)
{
#if NTRUGEN_LE && NTRUGEN_UNALIGNED
	*(uint16_t *)dst = x;
#else
	uint8_t *buf = dst;
	buf[0] = (uint8_t)x;
	buf[1] = (uint8_t)(x >> 8);
#endif
}

static inline uint32_t
dec32le(const void *src)
{
#if NTRUGEN_LE && NTRUGEN_UNALIGNED
	return *(const uint32_t *)src;
#else
	const uint8_t *buf = src;
	return (uint32_t)buf[0]
		| ((uint32_t)buf[1] << 8)
		| ((uint32_t)buf[2] << 16)
		| ((uint32_t)buf[3] << 24);
#endif
}

static inline void
enc32le(void *dst, uint32_t x)
{
#if NTRUGEN_LE && NTRUGEN_UNALIGNED
	*(uint32_t *)dst = x;
#else
	uint8_t *buf = dst;
	buf[0] = (uint8_t)x;
	buf[1] = (uint8_t)(x >> 8);
	buf[2] = (uint8_t)(x >> 16);
	buf[3] = (uint8_t)(x >> 24);
#endif
}

static inline uint64_t
dec64le(const void *src)
{
#if NTRUGEN_LE && NTRUGEN_UNALIGNED
	return *(const uint64_t *)src;
#else
	const uint8_t *buf = src;
	return (uint64_t)buf[0]
		| ((uint64_t)buf[1] << 8)
		| ((uint64_t)buf[2] << 16)
		| ((uint64_t)buf[3] << 24)
		| ((uint64_t)buf[4] << 32)
		| ((uint64_t)buf[5] << 40)
		| ((uint64_t)buf[6] << 48)
		| ((uint64_t)buf[7] << 56);
#endif
}

static inline void
enc64le(void *dst, uint64_t x)
{
#if NTRUGEN_LE && NTRUGEN_UNALIGNED
	*(uint64_t *)dst = x;
#else
	uint8_t *buf = dst;
	buf[0] = (uint8_t)x;
	buf[1] = (uint8_t)(x >> 8);
	buf[2] = (uint8_t)(x >> 16);
	buf[3] = (uint8_t)(x >> 24);
	buf[4] = (uint8_t)(x >> 32);
	buf[5] = (uint8_t)(x >> 40);
	buf[6] = (uint8_t)(x >> 48);
	buf[7] = (uint8_t)(x >> 56);
#endif
}

/* ==================================================================== */
/*
 * Modular arithmetics.
 *
 * We implement computations modulo some small integer p with the following
 * characteristics:
 *
 *   (4/3)*2^30 < p < 2^31       (this implies that 2*p < 2^32 < 3*p)
 *   p-1 is a multiple of 2048
 *
 * Operands are held in 32-bit values (uint32_t). We define R = 2^32 mod p.
 *
 * Values modulo p are 32-bit integers (uint32_t type) in the 0 to p-1 range.
 * Montgomery representation of an element x of Z_p is the value x*R mod p
 * (also in the 0 to p-1 range). Montgomery multiplication of x and y
 * computes x*y/R mod p (thus, Montgomery multiplication of the Montgomery
 * representations of x and y outputs the Montgomery representation of the
 * product x*y, since (x*R)*(y*R)/R = (x*y)*R). In general, values are not
 * kept in Montgomery representation, unless explicitly specified.
 *
 * The "signed normalized" value of x modulo p is the unique integer v
 * in the -(p-1)/2 to +(p-1)/2 range such that x = v mod p.
 *
 * The PRIMES[] array contains the largest of such primes, in
 * descending order. Six values are provided for each prime p:
 *    p     modulus
 *    p0i   -1/p mod 2^32
 *    R2    2^64 mod p
 *    g     a primitive 2048-th root of 1 modulo p (i.e. g^1024 = -1 mod p)
 *    ig    1/g mod p
 *    s     inverse mod p of the product of the previous primes
 * Values g, ig and s are in Montgomery representation.
 * R2 is used to convert values to Montgomery representations
 * (with montymul(x, R2) = x*R2/R = x*R mod p). g and ig are used to
 * generate the tables used for NTT and inverse NTT. Value s supports
 * reconstruction of a big integer in RNS representation with the CRT.
 * The product of all the primes in the PRIMES[] array is about 2^10012.25
 * and thus appropriate for big integers (in RNS) of up to 10000 bits.
 *
 * Polynomials over Z_p are considered modulo X^n+1 for n a power of two
 * between 2 and 1024 (inclusive). The degree n is provided as parameter
 * 'logn' in the 1 to 10 range (with n = 2^logn). The polynomial
 * coefficients are consecutive in RAM, in ascending degree order. The
 * NTT representation of such a polynomial is the evaluation of
 * the polynomial over the roots of X^n+1 (which are (g^(1024/n))^(2*i+1)
 * for integers i = 0 to n-1).
 *
 * Non-prime modulus
 * -----------------
 *
 * This code also works for some non-prime moduli p. If p is a product
 * p = p_1*p_2*...*p_k, such that each p_i is an odd prime, and p is in
 * the allowed range ((4/3)*2^30 to 2^31), then all operations work and
 * really compute things modulo each p_i simultaneously (through the
 * CRT). In particular, we can compute modulo q = 12289 by using (for
 * instance) p = 2013442049 = 163841*q (value 163841 is itself prime,
 * and 163840 is a multiple of 2048, which is not actually needed if we
 * are only interested in correct computations modulo q).
 */

/*
 * Expand the top bit of value x into a full 32-bit mask (i.e. return
 * 0xFFFFFFFF if x >= 0x80000000, or 0x00000000 otherwise).
 */
static inline uint32_t
tbmask(uint32_t x)
{
	return (uint32_t)(*(int32_t *)&x >> 31);
}

/*
 * Get v mod p in the 0 to p-1 range; input v must be in the -(p-1) to +(p-1)
 * range.
 */
static inline uint32_t
mp_set(int32_t v, uint32_t p)
{
	uint32_t w = (uint32_t)v;
	return w + (p & tbmask(w));
}

#if NTRUGEN_AVX2
TARGET_AVX2
static inline __m256i
mp_set_x8(__m256i yv, __m256i yp)
{
	return _mm256_add_epi32(yv, _mm256_and_si256(yp,
		_mm256_srai_epi32(yv, 31)));
}
#endif // NTRUGEN_AVX2

/*
 * Get the signed normalized value of x mod p.
 */
static inline int32_t
mp_norm(uint32_t x, uint32_t p)
{
	uint32_t w = x - (p & tbmask((p >> 1) - x));
	return *(int32_t *)&w;
}

#if NTRUGEN_AVX2
TARGET_AVX2
static inline __m256i
mp_norm_x8(__m256i yv, __m256i yp, __m256i yhp)
{
	return _mm256_sub_epi32(yv, _mm256_and_si256(yp,
		_mm256_cmpgt_epi32(yv, yhp)));
}
#endif // NTRUGEN_AVX2

#if 0 /* unused */
/*
 * Compute p0i = -1/p mod 2^32.
 */
static inline uint32_t
mp_ninv32(uint32_t p)
{
	uint32_t y = 2 - p;
	y *= 2 - p * y;
	y *= 2 - p * y;
	y *= 2 - p * y;
	y *= 2 - p * y;
	return -y;
}
#endif

/*
 * Compute R = 2^32 mod p.
 */
static inline uint32_t
mp_R(uint32_t p)
{
	/*
	 * Since 2*p < 2^32 < 3*p, we just subtract 2*p from 2^32.
	 */
	return -(p << 1);
}

/*
 * Compute R/2 = 2^31 mod p.
 */
static inline uint32_t
mp_hR(uint32_t p)
{
	/*
	 * Since p < 2^31 < (3/2)*p, we just subtract p from 2^31.
	 */
	return ((uint32_t)1 << 31) - p;
}

/*
 * Addition modulo p.
 */
static inline uint32_t
mp_add(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t d = a + b - p;
	return d + (p & tbmask(d));
}

#if NTRUGEN_AVX2
TARGET_AVX2
static inline __m256i
mp_add_x8(__m256i ya, __m256i yb, __m256i yp)
{
	__m256i yd = _mm256_sub_epi32(_mm256_add_epi32(ya, yb), yp);
	return _mm256_add_epi32(yd, _mm256_and_si256(yp,
		_mm256_srai_epi32(yd, 31)));
}
#endif // NTRUGEN_AVX2

/*
 * Subtraction modulo p.
 */
static inline uint32_t
mp_sub(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t d = a - b;
	return d + (p & tbmask(d));
}

#if NTRUGEN_AVX2
TARGET_AVX2
static inline __m256i
mp_sub_x8(__m256i ya, __m256i yb, __m256i yp)
{
	__m256i yd = _mm256_sub_epi32(ya, yb);
	return _mm256_add_epi32(yd, _mm256_and_si256(yp,
		_mm256_srai_epi32(yd, 31)));
}
#endif // NTRUGEN_AVX2

/*
 * Halving modulo p.
 */
static inline uint32_t
mp_half(uint32_t a, uint32_t p)
{
	return (a + (p & -(a & 1))) >> 1;
}

#if NTRUGEN_AVX2
TARGET_AVX2
static inline __m256i
mp_half_x8(__m256i ya, __m256i yp)
{
	return _mm256_srli_epi32(
		_mm256_add_epi32(ya, _mm256_and_si256(yp,
			_mm256_sub_epi32(_mm256_setzero_si256(),
			_mm256_and_si256(ya, _mm256_set1_epi32(1))))), 1);
}
#endif // NTRUGEN_AVX2

/*
 * Montgomery multiplication modulo p.
 *
 * Reduction computes (a*b + w*p)/(2^32) for some w <= 2^(32-1);
 * then p is conditionally subtracted. This process works as long as:
 *    (a*b + p*(2^32-1))/(2^32) <= 2*p-1
 * which holds if:
 *    a*b <= p*2^32 - 2^32 + p
 * This works if both a and b are proper integers modulo p (in the 0 to p-1
 * range), but also if, for instance, a is an integer modulo p, and b is an
 * arbitrary 32-bit integer.
 */
static inline uint32_t
mp_montymul(uint32_t a, uint32_t b, uint32_t p, uint32_t p0i)
{
	uint64_t z = (uint64_t)a * (uint64_t)b;
	uint32_t w = (uint32_t)z * p0i;
	uint32_t d = (uint32_t)((z + (uint64_t)w * (uint64_t)p) >> 32) - p;
	return d + (p & tbmask(d));
}

#if NTRUGEN_AVX2
/*
 * Input:
 *    ya = a0 : XX : a1 : XX : a2 : XX : a3 : XX
 *    yb = b0 : XX : b1 : XX : b2 : XX : b3 : XX
 * Output:
 *    mm(a0,b0) : 00 : mm(a1,b1) : 00 : mm(a2,b2) : 00 : mm(a3,b3) : 00
 */
TARGET_AVX2
static inline __m256i
mp_montymul_x4(__m256i ya, __m256i yb, __m256i yp, __m256i yp0i)
{
	__m256i yd = _mm256_mul_epu32(ya, yb);
	__m256i ye = _mm256_mul_epu32(yd, yp0i);
	ye = _mm256_mul_epu32(ye, yp);
	yd = _mm256_srli_epi64(_mm256_add_epi64(yd, ye), 32);
	yd = _mm256_sub_epi32(yd, yp);
	return _mm256_add_epi32(yd, _mm256_and_si256(yp,
		_mm256_srai_epi32(yd, 31)));
}

TARGET_AVX2
static inline __m256i
mp_montymul_x8(__m256i ya, __m256i yb, __m256i yp, __m256i yp0i)
{
	/* yd0 <- a0*b0 : a2*b2 (+high lane) */
	__m256i yd0 = _mm256_mul_epu32(ya, yb);
	/* yd1 <- a1*b1 : a3*b3 (+high lane) */
	__m256i yd1 = _mm256_mul_epu32(
		_mm256_srli_epi64(ya, 32),
		_mm256_srli_epi64(yb, 32));

	__m256i ye0 = _mm256_mul_epu32(yd0, yp0i);
	__m256i ye1 = _mm256_mul_epu32(yd1, yp0i);
	ye0 = _mm256_mul_epu32(ye0, yp);
	ye1 = _mm256_mul_epu32(ye1, yp);
	yd0 = _mm256_add_epi64(yd0, ye0);
	yd1 = _mm256_add_epi64(yd1, ye1);

	/* yf0 <- lo(d0) : lo(d1) : hi(d0) : hi(d1) (+high lane) */
	__m256i yf0 = _mm256_unpacklo_epi32(yd0, yd1);
	/* yf1 <- lo(d2) : lo(d3) : hi(d2) : hi(d3) (+high lane) */
	__m256i yf1 = _mm256_unpackhi_epi32(yd0, yd1);
	/* yg <- hi(d0) : hi(d1) : hi(d2) : hi(d3) (+high lane) */
	__m256i yg = _mm256_unpackhi_epi64(yf0, yf1);
	/*
	 * Alternate version (instead of the three unpack above) but it
	 * seems to be slightly slower.
	__m256i yg = _mm256_blend_epi32(_mm256_srli_epi64(yd0, 32), yd1, 0xAA);
	 */

	yg = _mm256_sub_epi32(yg, yp);
	return _mm256_add_epi32(yg, _mm256_and_si256(yp,
		_mm256_srai_epi32(yg, 31)));
}
#endif // NTRUGEN_AVX2

/*
 * Compute 2^(31*e) mod p.
 */
static inline uint32_t
mp_Rx31(unsigned e, uint32_t p, uint32_t p0i, uint32_t R2)
{
	/* x <- 2^63 mod p = Montgomery representation of 2^31 */
	uint32_t x = mp_half(R2, p);
	uint32_t d = 1;
	for (;;) {
		if ((e & 1) != 0) {
			d = mp_montymul(d, x, p, p0i);
		}
		e >>= 1;
		if (e == 0) {
			return d;
		}
		x = mp_montymul(x, x, p, p0i);
	}
}

/*
 * Division modulo p (x = dividend, y = divisor).
 * This code uses a constant-time binary GCD, which also works for a
 * non-prime modulus p (contrary to Fermat's Little Theorem). If the
 * divisor is not invertible modulo p, then 0 is returned.
 */
#define mp_div   Zn(mp_div)
uint32_t mp_div(uint32_t x, uint32_t y, uint32_t p);

#if NTRUGEN_AVX2
#define mp_div_x8   Zn(mp_div_x8)
TARGET_AVX2
__m256i mp_div_x8(__m256i ynum, __m256i yden, __m256i yp);
#endif // NTRUGEN_AVX2

/*
 * Compute the roots for NTT; given g (primitive 2048-th root of 1 modulo p),
 * this fills gm[] and igm[] with powers of g and 1/g:
 *    gm[rev(i)] = g^i mod p              (in Montgomery representation)
 *    igm[rev(i)] = (1/2)*(1/g)^i mod p   (in Montgomery representation)
 * rev() is the bit-reversal function over 10 bits. The arrays gm[] and igm[]
 * are filled only up to n = 2^logn values. Roots g and ig must be provided
 * in Montgomery representation.
 */
#define mp_mkgmigm   Zn(mp_mkgmigm)
void mp_mkgmigm(unsigned logn, uint32_t *restrict gm, uint32_t *restrict igm,
	uint32_t g, uint32_t ig, uint32_t p, uint32_t p0i);

/*
 * Like mp_mkgmigm(), but computing only gm[].
 */
#define mp_mkgm   Zn(mp_mkgm)
void mp_mkgm(unsigned logn, uint32_t *restrict gm,
	uint32_t g, uint32_t p, uint32_t p0i);

/*
 * A variant of mp_mkgm(), specialized for logn = 7, and g being a
 * 256-th root of 1, not a 2048-th root of 1.
 */
#define mp_mkgm7   Zn(mp_mkgm7)
void mp_mkgm7(uint32_t *restrict gm, uint32_t g, uint32_t p, uint32_t p0i);

/*
 * Like mp_mkgmigm(), but computing only igm[].
 */
#define mp_mkigm   Zn(mp_mkigm)
void mp_mkigm(unsigned logn, uint32_t *restrict igm,
	uint32_t ig, uint32_t p, uint32_t p0i);

/*
 * Compute the NTT over a polynomial. The polynomial a[] is modified in-place.
 */
#define mp_NTT   Zn(mp_NTT)
void mp_NTT(unsigned logn, uint32_t *restrict a, const uint32_t *restrict gm,
	uint32_t p, uint32_t p0i);

/*
 * Compute the inverse NTT over a polynomial. The polynomial a[] is modified
 * in-place.
 */
#define mp_iNTT   Zn(mp_iNTT)
void mp_iNTT(unsigned logn, uint32_t *restrict a, const uint32_t *restrict igm,
	uint32_t p, uint32_t p0i);

/*
 * Precomputed small primes. Enough values are provided to allow
 * computations in RNS representation over big integers up to 10000 bits.
 */
typedef struct {
	uint32_t p;
	uint32_t p0i;
	uint32_t R2;
	uint32_t g;
	uint32_t ig;
	uint32_t s;
} small_prime;
#define PRIMES   Zn(PRIMES)
extern const small_prime PRIMES[];

/* ==================================================================== */
/*
 * Custom bignum implementation.
 *
 * Big integers are represented as sequences of 32-bit integers; the
 * integer values are not necessarily consecutive in RAM (a dynamically
 * provided "stride" value is added to the current word pointer, to get
 * to the next word). The "len" parameter qualifies the number of words.
 *
 * Normal representation uses 31-bit limbs; each limb is stored in a
 * 32-bit word, with the top bit (31) always cleared. Limbs are in
 * low-to-high order. Signed integers use two's complement (hence, bit 30
 * of the last limb is the sign bit).
 *
 * RNS representation of a big integer x is the sequence of values
 * x modulo p, for the primes p defined in the PRIMES[] array.
 */

/*
 * Mutiply the provided big integer m with a small value x. The big
 * integer must have stride 1.
 * This function assumes that x < 2^31. The carry word is returned.
 */
#define zint_mul_small   Zn(zint_mul_small)
uint32_t zint_mul_small(uint32_t *m, size_t len, uint32_t x);

/*
 * Reduce a big integer d modulo a small integer p.
 * Rules:
 *  d is unsigned
 *  p is prime
 *  2^30 < p < 2^31
 *  p0i = -(1/p) mod 2^31
 *  R2 = 2^64 mod p
 */
#define zint_mod_small_unsigned   Zn(zint_mod_small_unsigned)
uint32_t zint_mod_small_unsigned(const uint32_t *d, size_t len, size_t stride,
	uint32_t p, uint32_t p0i, uint32_t R2);

#if NTRUGEN_AVX2
#define zint_mod_small_unsigned_x8 Zn(zint_mod_small_unsigned_x8)
TARGET_AVX2
__m256i zint_mod_small_unsigned_x8(
	const uint32_t *d, size_t len, size_t stride,
	__m256i yp, __m256i yp0i, __m256i yR2);
#endif // NTRUGEN_AVX2

/*
 * Similar to zint_mod_small_unsigned(), except that d may be signed.
 * Extra parameter is Rx = 2^(31*len) mod p.
 */
static inline uint32_t
zint_mod_small_signed(const uint32_t *d, size_t len, size_t stride,
	uint32_t p, uint32_t p0i, uint32_t R2, uint32_t Rx)
{
	if (len == 0) {
		return 0;
	}
	uint32_t z = zint_mod_small_unsigned(d, len, stride, p, p0i, R2);
	z = mp_sub(z, Rx & -(d[(len - 1) * stride] >> 30), p);
	return z;
}

#if NTRUGEN_AVX2
TARGET_AVX2
static inline __m256i
zint_mod_small_signed_x8(const uint32_t *d, size_t len, size_t stride,
	__m256i yp, __m256i yp0i, __m256i yR2, __m256i yRx)
{
	if (len == 0) {
		return _mm256_setzero_si256();
	}
	__m256i yz = zint_mod_small_unsigned_x8(d, len, stride, yp, yp0i, yR2);
	__m256i yl = _mm256_loadu_si256((__m256i *)(d + (len - 1) * stride));
	__m256i ym = _mm256_sub_epi32(_mm256_setzero_si256(),
		_mm256_srli_epi32(yl, 30));
	yz = mp_sub_x8(yz, _mm256_and_si256(yRx, ym), yp);
	return yz;
}
#endif // NTRUGEN_AVX2

/*
 * Add s*a to d. d and a initially have length 'len' words; the new d
 * has length 'len+1' words. 's' must fit on 31 bits. d[] and a[] must
 * not overlap. d uses stride dstride, while a has stride 1.
 */
#define zint_add_mul_small   Zn(zint_add_mul_small)
void zint_add_mul_small(uint32_t *restrict d, size_t len, size_t dstride,
	const uint32_t *restrict a, uint32_t s);

#if NTRUGEN_AVX2
/*
 * Like zint_add_mul_small(), except that it handles eight integers in
 * parallel:
 *    d0 <- d0 + s0*a
 *    d1 <- d1 + s1*a
 *     ...
 *    d7 <- d7 + s7*a
 */
#define zint_add_mul_small_x8   Zn(zint_add_mul_small_x8)
TARGET_AVX2
void zint_add_mul_small_x8(uint32_t *restrict d, size_t len, size_t dstride,
	const uint32_t *restrict a, __m256i ys);
#endif // NTRUGEN_AVX2

/*
 * Normalize a modular integer around 0: if x > p/2, then x is replaced
 * with x - p (signed encoding with two's complement); otherwise, x is
 * untouched. The two integers x and p are encoded over the same length;
 * x has stride xstride, while p has stride 1.
 */
#define zint_norm_zero   Zn(zint_norm_zero)
void zint_norm_zero(uint32_t *restrict x, size_t len, size_t xstride,
	const uint32_t *restrict p);

/*
 * Rebuild integers from their RNS representation. There are 'num_sets' sets
 * of 'n' integers. Within each set, the n integers are interleaved,
 * so that words of a given integer occur every n slots in RAM (i.e. each
 * integer has stride 'n'). The sets are consecutive in RAM.
 *
 * If "normalize_signed" is non-zero, then the output values are
 * normalized to the -m/2..m/2 interval (where m is the product of all
 * small prime moduli); two's complement is used for negative values.
 * If "normalize_signed" is zero, then the output values are all
 * in the 0..m-1 range.
 *
 * tmp[] must have room for xlen words.
 */
#define zint_rebuild_CRT   Zn(rebuild_CRT)
void zint_rebuild_CRT(uint32_t *restrict xx, size_t xlen, size_t n,
	size_t num_sets, int normalize_signed, uint32_t *restrict tmp);

/*
 * Negate a big integer conditionally: value a is replaced with -a if
 * and only if ctl = 1. Control value ctl must be 0 or 1. The integer
 * has stride 1.
 */
#define zint_negate   Zn(zint_negate)
void zint_negate(uint32_t *a, size_t len, uint32_t ctl);

/*
 * Get the number of leading zeros in a 32-bit value.
 */
TARGET_AVX2
static inline unsigned
lzcnt(uint32_t x)
{
#if NTRUGEN_AVX2
	/*
	 * All AVX2-capable CPUs have lzcnt.
	 */
	return _lzcnt_u32(x);
#else // NTRUGEN_AVX2
	uint32_t m = tbmask((x >> 16) - 1);
	uint32_t s = m & 16;
	x = (x >> 16) ^ (m & (x ^ (x >> 16)));
	m = tbmask((x >>  8) - 1);
	s |= m &  8;
	x = (x >>  8) ^ (m & (x ^ (x >>  8)));
	m = tbmask((x >>  4) - 1);
	s |= m &  4;
	x = (x >>  4) ^ (m & (x ^ (x >>  4)));
	m = tbmask((x >>  2) - 1);
	s |= m &  2;
	x = (x >>  2) ^ (m & (x ^ (x >>  2)));

	/*
	 * At this point, x fits on 2 bits. Number of leading zeros is
	 * then:
	 *    x = 0   -> 2
	 *    x = 1   -> 1
	 *    x = 2   -> 0
	 *    x = 3   -> 0
	 */
	return (unsigned)(s + ((2 - x) & tbmask(x - 3)));
#endif // NTRUGEN_AVX2
}

/*
 * Identical to lzcnt(), except that the caller makes sure that the
 * operand is non-zero. On (old-ish) x86 systems, this function could be
 * specialized with the bsr opcode (which does not support a zero input).
 */
#define lzcnt_nonzero   lzcnt

/*
 * Compute a GCD between two positive big integers x and y. The two
 * integers must be odd. Returned value is 1 if the GCD is 1, 0
 * otherwise. When 1 is returned, arrays u and v are filled with values
 * such that:
 *   0 <= u <= y
 *   0 <= v <= x
 *   x*u - y*v = 1
 * x[] and y[] are unmodified. Both input values must have the same
 * encoded length. Temporary array must be large enough to accommodate 4
 * extra values of that length. Arrays u, v and tmp may not overlap with
 * each other, or with either x or y. All integers use stride 1.
 */
#define zint_bezout   Zn(zint_bezout)
int zint_bezout(uint32_t *restrict u, uint32_t *restrict v,
	const uint32_t *restrict x, const uint32_t *restrict y,
	size_t len, uint32_t *restrict tmp);

/*
 * Add k*(2^sc)*y to x. The result is assumed to fit in the array of
 * size xlen (truncation is applied if necessary).
 * Scale factor sc is provided as sch and scl, such that:
 *    sch = sc / 31
 *    scl = sc % 31  (in the 0..30 range)
 * xlen MUST NOT be lower than ylen; however, it is allowed that
 * xlen is greater than ylen.
 *
 * x[] and y[] are both signed integers, using two's complement for
 * negative values. They both use the same stride ('stride' parameter).
 */
#define zint_add_scaled_mul_small   Zn(zint_add_scaled_mul_small)
void zint_add_scaled_mul_small(uint32_t *restrict x, size_t xlen,
	const uint32_t *restrict y, size_t ylen, size_t stride,
	int32_t k, uint32_t sch, uint32_t scl);

/*
 * Subtract y*2^sc from x. This is a specialized version of
 * zint_add_scaled_mul_small(), with multiplier k = -1.
 */
#define zint_sub_scaled   Zn(zint_sub_scaled)
void zint_sub_scaled(uint32_t *restrict x, size_t xlen,
	const uint32_t *restrict y, size_t ylen, size_t stride,
	uint32_t sch, uint32_t scl);

/* ====================================================================== */
/*
 * Fixed-point numbers.
 *
 * For FFT and other computations with approximations, we use a fixed-point
 * format over 64 bits; the top 32 bits are the integral part, and the low
 * 32 bits are the fractional part.
 */

/*
 * We wrap the type into a struct in order to detect any attempt at using
 * arithmetic operators on values directly. Since all functions are inline,
 * the compiler will be able to remove the wrapper, which will then have
 * no runtime cost.
 */
typedef struct {
	uint64_t v;
} fxr;

#define FXR(x)   { (x) }

static inline fxr
fxr_of(int32_t j)
{
	fxr x;

	x.v = (uint64_t)j << 32;
	return x;
}

static inline fxr
fxr_of_scaled32(uint64_t t)
{
	fxr x;

	x.v = t;
	return x;
}

static inline fxr
fxr_add(fxr x, fxr y)
{
	x.v += y.v;
	return x;
}

static inline fxr
fxr_sub(fxr x, fxr y)
{
	x.v -= y.v;
	return x;
}

static inline fxr
fxr_double(fxr x)
{
	x.v <<= 1;
	return x;
}

static inline fxr
fxr_neg(fxr x)
{
	x.v = -x.v;
	return x;
}

static inline fxr
fxr_abs(fxr x)
{
	x.v -= (x.v << 1) & (uint64_t)(*(int64_t *)&x.v >> 63);
	return x;
}

static inline fxr
fxr_mul(fxr x, fxr y)
{
#if defined __GNUC__ && defined __SIZEOF_INT128__
	__int128 z;

	z = (__int128)*(int64_t *)&x.v * (__int128)*(int64_t *)&y.v;
	x.v = (uint64_t)(z >> 32);
	return x;
#else
	int32_t xh, yh;
	uint32_t xl, yl;
	uint64_t z0, z1, z2, z3;

	xl = (uint32_t)x.v;
	yl = (uint32_t)y.v;
	xh = (int32_t)(*(int64_t *)&x.v >> 32);
	yh = (int32_t)(*(int64_t *)&y.v >> 32);
	z0 = ((uint64_t)xl * (uint64_t)yl) >> 32;
	z1 = (uint64_t)((int64_t)xl * (int64_t)yh);
	z2 = (uint64_t)((int64_t)yl * (int64_t)xh);
	z3 = (uint64_t)((int64_t)xh * (int64_t)yh) << 32;
	x.v = z0 + z1 + z2 + z3;
	return x;
#endif
}

#if NTRUGEN_AVX2
TARGET_AVX2
static inline __m256i
fxr_mul_x4(__m256i ya, __m256i yb)
{
	__m256i ya_hi = _mm256_srli_epi64(ya, 32);
	__m256i yb_hi = _mm256_srli_epi64(yb, 32);
	__m256i y1 = _mm256_mul_epu32(ya, yb);
	__m256i y2 = _mm256_mul_epu32(ya, yb_hi);
	__m256i y3 = _mm256_mul_epu32(ya_hi, yb);
	__m256i y4 = _mm256_mul_epu32(ya_hi, yb_hi);
	y1 = _mm256_srli_epi64(y1, 32);
	y4 = _mm256_slli_epi64(y4, 32);
	__m256i y5 = _mm256_add_epi64(
		_mm256_add_epi64(y1, y2),
		_mm256_add_epi64(y3, y4));
	__m256i yna = _mm256_srai_epi32(ya, 31);
	__m256i ynb = _mm256_srai_epi32(yb, 31);
	return _mm256_sub_epi64(y5,
		_mm256_add_epi64(
			_mm256_and_si256(_mm256_slli_epi64(yb, 32), yna),
			_mm256_and_si256(_mm256_slli_epi64(ya, 32), ynb)));
}
#endif // NTRUGEN_AVX2

static inline fxr
fxr_sqr(fxr x)
{
#if defined __GNUC__ && defined __SIZEOF_INT128__
	int64_t t;
	__int128 z;

	t = *(int64_t *)&x.v;
	z = (__int128)t * (__int128)t;
	x.v = (uint64_t)(z >> 32);
	return x;
#else
	int32_t xh;
	uint32_t xl;
	uint64_t z0, z1, z3;

	xl = (uint32_t)x.v;
	xh = (int32_t)(*(int64_t *)&x.v >> 32);
	z0 = ((uint64_t)xl * (uint64_t)xl) >> 32;
	z1 = (uint64_t)((int64_t)xl * (int64_t)xh);
	z3 = (uint64_t)((int64_t)xh * (int64_t)xh) << 32;
	x.v = z0 + (z1 << 1) + z3;
	return x;
#endif
}

#if NTRUGEN_AVX2
TARGET_AVX2
static inline __m256i
fxr_sqr_x4(__m256i ya)
{
	__m256i ya_hi = _mm256_srli_epi64(ya, 32);
	__m256i y1 = _mm256_mul_epu32(ya, ya);
	__m256i y2 = _mm256_mul_epu32(ya, ya_hi);
	__m256i y3 = _mm256_mul_epu32(ya_hi, ya_hi);
	y1 = _mm256_srli_epi64(y1, 32);
	y2 = _mm256_add_epi64(y2, y2);
	y3 = _mm256_slli_epi64(y3, 32);
	__m256i y4 = _mm256_add_epi64(_mm256_add_epi64(y1, y2), y3);
	return _mm256_sub_epi64(y4,
		_mm256_and_si256(_mm256_slli_epi64(ya, 33),
		_mm256_srai_epi32(ya, 31)));
}
#endif // NTRUGEN_AVX2

static inline int32_t
fxr_round(fxr x)
{
	x.v += 0x80000000ul;
	return (int32_t)(*(int64_t *)&x.v >> 32);
}

static inline fxr
fxr_div2e(fxr x, unsigned n)
{
	x.v += (((uint64_t)1 << n) >> 1);
	x.v = (uint64_t)(*(int64_t *)&x.v >> n);
	return x;
}

#if NTRUGEN_AVX2
TARGET_AVX2
static inline __m256i
fxr_half_x4(__m256i ya)
{
	const __m256i y1 = _mm256_set1_epi64x(1);
	const __m256i yh = _mm256_set1_epi64x((uint64_t)1 << 63);
	ya = _mm256_add_epi64(ya, y1);
	return _mm256_or_si256(
		_mm256_srli_epi64(ya, 1),
		_mm256_and_si256(ya, yh));
}
#endif // NTRUGEN_AVX2

static inline fxr
fxr_mul2e(fxr x, unsigned n)
{
	x.v <<= n;
	return x;
}

#define inner_fxr_div   Zn(inner_fxr_div)
uint64_t inner_fxr_div(uint64_t x, uint64_t y);

static inline fxr
fxr_inv(fxr x)
{
	x.v = inner_fxr_div((uint64_t)1 << 32, x.v);
	return x;
}

static inline fxr
fxr_div(fxr x, fxr y)
{
	x.v = inner_fxr_div(x.v, y.v);
	return x;
}

#if NTRUGEN_AVX2
#define fxr_div_x4   Zn(fxr_div_x4)
TARGET_AVX2 __m256i fxr_div_x4(__m256i yn, __m256i yd);

/*
 * Divide four values (n0..n3) by the same divisor (d).
 */
TARGET_AVX2
static inline void
fxr_div_x4_1(fxr *n0, fxr *n1, fxr *n2, fxr *n3, fxr d)
{
	__m256i yn = _mm256_setr_epi64x(n0->v, n1->v, n2->v, n3->v);
	__m256i yd = _mm256_set1_epi64x(d.v);
	union {
		__m256i y;
		uint64_t q[4];
	} z;
	z.y = fxr_div_x4(yn, yd);
	n0->v = z.q[0];
	n1->v = z.q[1];
	n2->v = z.q[2];
	n3->v = z.q[3];
}
#endif // NTRUGEN_AVX2

static inline int
fxr_lt(fxr x, fxr y)
{
	return *(int64_t *)&x.v < *(int64_t *)&y.v;
}

static const fxr fxr_zero = { 0 };
static const fxr fxr_sqrt2 = { 6074001000ull };

/*
 * A complex value.
 */
typedef struct {
	fxr re, im;
} fxc;

#define FXC(re, im)   { FXR(re), FXR(im) }

static inline fxc
fxc_add(fxc x, fxc y)
{
	x.re = fxr_add(x.re, y.re);
	x.im = fxr_add(x.im, y.im);
	return x;
}

static inline fxc
fxc_sub(fxc x, fxc y)
{
	x.re = fxr_sub(x.re, y.re);
	x.im = fxr_sub(x.im, y.im);
	return x;
}

static inline fxc
fxc_half(fxc x)
{
	x.re = fxr_div2e(x.re, 1);
	x.im = fxr_div2e(x.im, 1);
	return x;
}

static inline fxc
fxc_mul(fxc x, fxc y)
{
	/*
	 * We are computing r = (a + i*b)*(c + i*d) with:
	 *   z0 = a*c
	 *   z1 = b*d
	 *   z2 = (a + b)*(c + d)
	 *   r = (z0 - z1) + i*(z2 - (z0 + z1))
	 * Since the intermediate values are truncated to our precision,
	 * the imaginary value of r _may_ be slightly different from
	 * a*d + b*c (if we had calculated it directly). For full
	 * reproducibility, all implementations should use the formulas
	 * above.
	 */
	fxr z0 = fxr_mul(x.re, y.re);
	fxr z1 = fxr_mul(x.im, y.im);
	fxr z2 = fxr_mul(fxr_add(x.re, x.im), fxr_add(y.re, y.im));
	fxc z;
	z.re = fxr_sub(z0, z1);
	z.im = fxr_sub(z2, fxr_add(z0, z1));
	return z;
}

#if NTRUGEN_AVX2
TARGET_AVX2
static inline void
fxc_mul_x4(__m256i *yd_re, __m256i *yd_im,
	__m256i ya_re, __m256i ya_im, __m256i yb_re, __m256i yb_im)
{
	__m256i y0 = fxr_mul_x4(ya_re, yb_re);
	__m256i y1 = fxr_mul_x4(ya_im, yb_im);
	__m256i y2 = fxr_mul_x4(
		_mm256_add_epi64(ya_re, ya_im),
		_mm256_add_epi64(yb_re, yb_im));
	*yd_re = _mm256_sub_epi64(y0, y1);
	*yd_im = _mm256_sub_epi64(y2, _mm256_add_epi64(y0, y1));
}
#endif // NTRUGEN_AVX2

static inline fxc
fxc_conj(fxc x)
{
	x.im = fxr_neg(x.im);
	return x;
}

/*
 * In FFT representation, we keep only half of the coefficients, because
 * all our vectors are real in non-FFT representation; thus, the FFT
 * representation is redundant. For 0 <= k < n/2, f[k] contains the
 * real part of FFT coefficient k, and f[k + n/2] contains the imaginary
 * part of FFT coefficient k.
 */

/*
 * Convert a (real) vector to its FFT representation.
 */
#define vect_FFT   Zn(vect_FFT)
void vect_FFT(unsigned logn, fxr *f);

/*
 * Convert back from FFT representation into a real vector.
 */
#define vect_iFFT   Zn(vect_iFFT)
void vect_iFFT(unsigned logn, fxr *f);

/*
 * Set a vector d to the value of the small polynomial f.
 */
#define vect_set   Zn(vect_set)
void vect_set(unsigned logn, fxr *d, const int8_t *f);

/*
 * Add vector b to vector a. This works in both real and FFT representations.
 * Vectors a and b MUST NOT overlap.
 */
#define vect_add   Zn(vect_add)
void vect_add(unsigned logn, fxr *restrict a, const fxr *restrict b);

/*
 * Multiply vector a by the real constant c. This works in both real
 * and FFT representations.
 */
#define vect_mul_realconst   Zn(vect_mul_realconst)
void vect_mul_realconst(unsigned logn, fxr *a, fxr c);

/*
 * Multiply a vector by 2^e. This works in both real and FFT representations.
 */
#define vect_mul2e   Zn(vect_mul2e)
void vect_mul2e(unsigned logn, fxr *a, unsigned e);

/*
 * Multiply vector a by vector b. The vectors must be in FFT representation.
 * Vectors a and b MUST NOT overlap.
 */
#define vect_mul_fft   Zn(vect_mul_fft)
void vect_mul_fft(unsigned logn, fxr *restrict a, const fxr *restrict b);

/*
 * Convert a vector into its adjoint (in FFT representation).
 */
#define vect_adj_fft   Zn(vect_adj_fft)
void vect_adj_fft(unsigned logn, fxr *a);

/*
 * Multiply vector a by auto-adjoint vector b. The vectors must be in FFT
 * representation. Since the FFT representation of an auto-adjoint vector
 * contains only real number, the second half of b contains only zeros and
 * is not accessed by this function. Vectors a and b MUST NOT overlap.
 */
#define vect_mul_autoadj_fft   Zn(vect_mul_autoadj_fft)
void vect_mul_autoadj_fft(unsigned logn,
	fxr *restrict a, const fxr *restrict b);

/*
 * Divide vector a by auto-adjoint vector b. The vectors must be in FFT
 * representation. Since the FFT representation of an auto-adjoint vector
 * contains only real number, the second half of b contains only zeros and
 * is not accessed by this function. Vectors a and b MUST NOT overlap.
 */
#define vect_div_autoadj_fft   Zn(vect_div_autoadj_fft)
void vect_div_autoadj_fft(unsigned logn,
	fxr *restrict a, const fxr *restrict b);

/*
 * Compute d = (2^e)/(a*adj(a) + b*adj(b)). Polynomials are in FFT
 * representation. Since d is auto-adjoint, only its first half is set; the
 * second half is _implicitly_ zero (this function does not access the
 * second half of d). Vectors a, b and d MUST NOT overlap.
 */
#define vect_invnorm_fft   Zn(vect_invnorm_fft)
void vect_invnorm_fft(unsigned logn, fxr *restrict d,
	const fxr *restrict a, const fxr *restrict b, unsigned e);

/*
 * Compute d = (2^e)*adj(a)/(a*adj(a)), written back into a[]. Polynomials
 * are in FFT representation.
 */
#define vect_inv_mul2e_fft   Zn(vect_inv_mul2e_fft)
void vect_inv_mul2e_fft(unsigned logn, fxr *a, unsigned e);

/* ==================================================================== */
/*
 * Code for polynomials with integer coefficients.
 *
 * Polynomials use an interleaved in-memory representation:
 *
 *   There are n = 2^logn coefficients (degrees 0 to n-1).
 *   Each coefficient contains 'len' words (zint limbs or RNS).
 *   Each coefficient has a stride of 'n'.
 *   The first (lowest) words of the respective coefficients are consecutive.
 */

/*
 * Load a one-byte polynomial with reduction modulo p.
 */
#define poly_mp_set_small   Zn(poly_mp_set_small)
void poly_mp_set_small(unsigned logn, uint32_t *restrict d,
	const int8_t *restrict f, uint32_t p);

/*
 * Convert a polynomial in one-word normal representation (signed) into RNS
 * modulo the single prime p.
 */
#define poly_mp_set   Zn(poly_mp_set)
void poly_mp_set(unsigned logn, uint32_t *f, uint32_t p);

/*
 * Convert a polynomial in RNS (modulo a single prime p) into one-word
 * normal representation (signed).
 */
#define poly_mp_norm   Zn(poly_mp_norm)
void poly_mp_norm(unsigned logn, uint32_t *f, uint32_t p);

/*
 * Convert a polynomial to small integers. Source values are supposed
 * to be normalized (signed). Returned value is 0 if any of the
 * coefficients exceeds the provided limit (in absolute value); on
 * success, 1 is returned.
 *
 * In case of failure, the function returns earlier; this does not
 * break constant-time discipline as long as a failure implies that the
 * (f,g) polynomials are discarded.
 */
#define poly_big_to_small   Zn(poly_big_to_small)
int poly_big_to_small(unsigned logn, int8_t *restrict d,
	const uint32_t *restrict s, int lim);

/*
 * Get the maximum bit length of all coefficients of a polynomial. Each
 * coefficient has size flen words.
 *
 * The bit length of a big integer is defined to be the length of the
 * minimal binary representation, using two's complement for negative
 * values, and excluding the sign bit. This definition implies that
 * if x = 2^k, then x has bit length k but -x has bit length k-1. For
 * non powers of two, x and -x have the same bit length.
 *
 * This function is constant-time with regard to coefficient values and
 * the returned bit length.
 */
#define poly_max_bitlength   Zn(poly_max_bitlength)
uint32_t poly_max_bitlength(unsigned logn, const uint32_t *f, size_t flen);

/*
 * Compute q = x / 31 and r = x % 31 for an unsigned integer x. This
 * macro is constant-time and works for values x up to 63487 (inclusive).
 */
#define DIVREM31(q, r, x)  { \
		uint32_t divrem31_q, divrem31_x; \
		divrem31_x = (x); \
		divrem31_q = (uint32_t)(divrem31_x * (uint32_t)67651) >> 21; \
		(q) = divrem31_q; \
		(r) = divrem31_x - 31 * divrem31_q; \
	} while (0)

/*
 * Convert a polynomial to fixed-point approximations, with scaling.
 * For each coefficient x, the computed approximation is x/2^sc.
 * This function assumes that |x| < 2^(30+sc). The length of each
 * coefficient must be less than 2^24 words.
 *
 * This function is constant-time with regard to the coefficient values
 * and to the scaling factor.
 */
#define poly_big_to_fixed   Zn(poly_big_to_fixed)
void poly_big_to_fixed(unsigned logn, fxr *restrict d,
	const uint32_t *restrict f, size_t len, uint32_t sc);

/*
 * Subtract k*f from F, where F, f and k are polynomials modulo X^n+1.
 * Coefficients of polynomial k are small integers (signed values in the
 * -2^31..+2^31 range) scaled by 2^sc.
 *
 * This function implements the basic quadratic multiplication algorithm,
 * which is efficient in space (no extra buffer needed) but slow at
 * high degree.
 */
#define poly_sub_scaled   Zn(poly_sub_scaled)
void poly_sub_scaled(unsigned logn,
	uint32_t *restrict F, size_t Flen,
	const uint32_t *restrict f, size_t flen,
	const int32_t *restrict k, uint32_t sc);

/*
 * Subtract k*f from F. Coefficients of polynomial k are small integers
 * (signed values in the -2^31..+2^31 range) scaled by 2^sc. Polynomial f
 * MUST be in RNS+NTT over flen+1 words (even though f itself would fit on
 * flen words); polynomial F MUST be in plain representation.
 */
#define poly_sub_scaled_ntt   Zn(poly_sub_scaled_ntt)
void
poly_sub_scaled_ntt(unsigned logn, uint32_t *restrict F, size_t Flen,
	const uint32_t *restrict f, size_t flen,
	const int32_t *restrict k, uint32_t sc, uint32_t *restrict tmp);

/*
 * depth = 1
 * logn = logn_top - depth
 * Inputs:
 *    F       polynomial of degree 2^logn, plain integer representation (FGlen)
 *    FGlen   size of each coefficient of F (must be 1 or 2)
 *    f       polynomial of degree 2^logn_top, small coefficients
 *    k       polynomial of degree 2^logn (plain, 32-bit)
 *    sc      scaling logarithm (public value)
 *    tmp     temporary with room at least max(FGlen, 2^logn_top) words
 * Operation:
 *    F <- F - (2^sc)*k*ft
 * with ft being the degree-n polynomial corresponding to f
 * It is assumed that the result fits.
 *
 * WARNING: polynomial k is consumed in the process.
 *
 * This function uses 3*n words in tmp[].
 */
#define poly_sub_kf_scaled_depth1   Zn(poly_sub_kf_scaled_depth1)
void poly_sub_kf_scaled_depth1(unsigned logn_top,
	uint32_t *restrict F, size_t FGlen,
	uint32_t *restrict k, uint32_t sc,
	const int8_t *restrict f,
	uint32_t *restrict tmp);

/*
 * Check whether the provided polynomial is invertible modulo X^n+1
 * and modulo some small prime r which is such that r-1 is a multiple
 * of 2048. Parameters:
 *    logn     degree logarithm
 *    f        polynomial to test
 *    r        small prime to use
 *    p        p = r*t, such that t is prime and (4/3)*2^30 < p < 2^31
 *    p0i      -1/p mod 2^32
 *    s        s'*2^32 mod p, for some s' such that s'^1024 = -1 mod r
 *    rm, rs   division by r parameters
 *    tmp      temporary buffer
 * rm and rs must be such that floor(x*rm/(2^rs)) == floor(x/r) for all
 * x in [0..p-1]; such values always exist.
 *
 * Return value: 1 if f is invertible, 0 otherwise.
 *
 * RAM USAGE: 2*n words
 */
#define poly_is_invertible   Zn(poly_is_invertible)
int
poly_is_invertible(unsigned logn, const int8_t *restrict f,
	uint32_t p, uint32_t p0i, uint32_t s,
	uint32_t r, uint32_t rm, unsigned rs, uint32_t *restrict tmp);

/*
 * Similar to poly_is_invertible(), except that we test two prime moduli
 * at the same time. We have p = r1*r2*t, with both r1-1 and r2-1 being
 * multiples of 2048, and s is a 2048-th root of 1 modulo both r1 and r2
 * (s is in Montgomery representation).
 */
#define poly_is_invertible_ext   Zn(poly_is_invertible_ext)
int
poly_is_invertible_ext(unsigned logn, const int8_t *restrict f,
	uint32_t r1, uint32_t r2, uint32_t p, uint32_t p0i, uint32_t s,
	uint32_t r1m, unsigned r1s, uint32_t r2m, unsigned r2s,
	uint32_t *restrict tmp);

/* ==================================================================== */
/*
 * For a Gaussian distribution centred on zero and with standard deviation
 * sigma, we define:
 *    D(x) = exp(-(x^2)/(2*sigma^2))
 * The probability of integer k to be selected is then D(k)/(sum_j D(j)).
 *
 * We then define probabilities:
 *    P(k) = (\sum_{j<=k} D(j))/(\sum_{all j} D(j))
 * i.e. P(k) is the probability of the selected integer to be at most k.
 * We scale and round these values into:
 *    Q(k) = round((2^16-1)*P(k))
 * (Note: 2^16-1, not 2^16, so that the max value is 65535)
 * We then store these values Q(k). The head values are all equal to 0,
 * the tail values are all equal to 65535; each table contains the values
 * Q(k) for k = -kmax to +(kmax-1), for the integer kmax such that
 * Q(-kmax) != 0 but Q(-kmax-1) = 0.
 *
 * Actual selection consists in the following:
 *    generate a random 16-bit integer x (with 0 <= x <= 65535)
 *    k = -kmax
 *    while (k < +kmax) and (Q(k) < x):
 *        k = k + 1
 *    return k
 * This process returns a value in the [-kmax; +kmax] range.
 *
 * The following deviations (sigma) are used, for various algorithms:
 *
 *    name          degree n       sigma
 *    Falcon-n      2 to 1024      1.17*sqrt(12289/n)
 *    BAT-128-256      256         sqrt(sqrt(2)*128/512)
 *    BAT-257-512      512         sqrt(sqrt(2)*257/1024)
 *    BAT-769-1024    1024         sqrt(sqrt(5*log(5)/24)*769/1024)
 *    Hawk-256         256         1.1
 *    Hawk-512         512         1.5
 *    Hawk-1024       1024         2.0
 *
 * All algorithms target f*G - g*F = q for some integer q; Falcon uses
 * q = 12289; BAT uses q = 128, 257 or 769 (depending on the degree n);
 * Hawk uses q = 1.
 *
 * Table format:
 * -------------
 *
 * The in-memory format for each table is the value of kmax, followed
 * by the 2*kmax values of Q(k) for k = -kmax to +kmax-1.
 *
 * Special case for Falcon reduced versions:
 * -----------------------------------------
 * Standard Falcon uses n = 512 or 1024. We include a table for n = 256
 * ("small" version, maybe useful for some use cases with very small
 * embedded systems?). We also support reduced versions (useful for
 * tests) for degrees n = 2 to 128, by using the n = 256 table
 * repeatedly. Namely, we use the table 256/n times, and add the values
 * together. This still follows the proper distribution, without
 * requiring any extra table. We have to enforce an arbitrary limit on
 * values, so that the resulting f and g are encodable in a way
 * compatible with the reference Falcon code; namely, for n <= 32, we
 * reject values which are not in [-127..+127]. This conditional
 * rejection means that the source RNG may have to be invoked a few more
 * times to get enough random material. This does not happen for the
 * degrees 64 to 1024, where the maximum output value (kmax) always fits
 * within the encoding format limit.
 */

/* BAT, q = 128, n = 256 -> kmax = 2 */
#define gauss_BAT_128_256   Zn(gauss_BAT_128_256)
extern const uint16_t gauss_BAT_128_256[];

/* BAT, q = 257, n = 512 -> kmax = 2 */
#define gauss_BAT_257_512   Zn(gauss_BAT_257_512)
extern const uint16_t gauss_BAT_257_512[];

/* BAT, q = 769, n = 1024 -> kmax = 3 */
#define gauss_BAT_769_1024   Zn(gauss_BAT_769_1024)
extern const uint16_t gauss_BAT_769_1024[];

/* Falcon, q = 12289, n = 256 -> kmax = 24 */
#define gauss_Falcon_256   Zn(gauss_Falcon_256)
extern const uint16_t gauss_Falcon_256[];

/* Falcon, q = 12289, n = 512 -> kmax = 17 */
#define gauss_Falcon_512   Zn(gauss_Falcon_512)
extern const uint16_t gauss_Falcon_512[];

/* Falcon, q = 12289, n = 1024 -> kmax = 12 */
#define gauss_Falcon_1024   Zn(gauss_Falcon_1024)
extern const uint16_t gauss_Falcon_1024[];

/*
 * For a given random source and selection table, generate polynomial f
 * This function ensures that the returned f has odd parity; if the
 * produced polynomial has even parity, it is discarded and the process
 * loops.
 */
#define gauss_sample_poly   Zn(gauss_sample_poly)
void gauss_sample_poly(unsigned logn, int8_t *f,
	const uint16_t *tab, ntrugen_rng rng, void *rng_context);

/*
 * Special sampling function for reduced versions: provided table should
 * be for degree n=256, and 256/(2^n) samples are added together for each
 * of the 2^n output values (re-sampling is performed when a value does
 * not fit in [-127..+127], so the amount of consumed randomness may vary;
 * minimum usage is 512 bytes).
 */
#define gauss_sample_poly_reduced   Zn(gauss_sample_poly_reduced)
void gauss_sample_poly_reduced(unsigned logn, int8_t *f,
	const uint16_t *tab, ntrugen_rng rng, void *rng_context);

/* see ng_inner.h */
#define poly_sqnorm   Zn(poly_sqnorm)
uint32_t poly_sqnorm(unsigned logn, const int8_t *f);

/* ==================================================================== */

/*
 * Lengths of values (big integers) in 31-bit limbs:
 *   max_bl_small[d]      max length of coeffs of (f,g) at depth d
 *   max_bl_small[logn]   max length of Res(f,phi) and Res(g,phi)
 *   max_bl_large[d]      max length of coeffs of unreduced (F,G) at depth d
 *   word_win[d]          number of top limbs to consider in (f,g) (Babai's NP)
 * Rules:
 *   max_bl_small[0] = 1
 *   max_bl_large[d] >= max_bl_small[d + 1]
 *   1 <= word_win[d] <= max_bl_small[d]
 * Additional rules to use the optimized depth0 function:
 *   max_bl_large[0] = 1
 *   max_bl_small[1] = 1
 *
 * q                target integer (f*G - g*F = q)
 * min_logn         minimum logn supported
 * max_logn         maximum logn supported
 * reduce_bits      assumed reduction per iteration (in bits)
 * coeff_FG_limit   maximum allowed value for coefficients of F and G
 * min_save_fg      minimum depth at which (f,g) can be saved temporarily
 */
typedef struct {
	uint32_t q;
	unsigned min_logn, max_logn;
	uint16_t max_bl_small[11];
	uint16_t max_bl_large[10];
	uint16_t word_win[10];
	uint32_t reduce_bits;
	uint8_t coeff_FG_limit[11];
	uint16_t min_save_fg[11];
} ntru_profile;

/* Error code: no error (so far) */
#define SOLVE_OK           0

/* Error code: GCD(Res(f,X^n+1), Res(g,X^n+1)) != 1 */
#define SOLVE_ERR_GCD      -1

/* Error code: reduction error (NTRU equation no longer fulfilled) */
#define SOLVE_ERR_REDUCE   -2

/* Error code: output (F,G) coefficients are off-limits */
#define SOLVE_ERR_LIMIT    -3

/*
 * Solve the NTRU equation for the provided (f,g).
 * The (F,G) solution (if found) is returned at the start of the tmp[]
 * array, as two consecutive int8_t[] values. Returned value is
 * SOLVE_OK (0) on success, a negative error code on failure.
 *
 * Note: if f is not invertible modulo X^n+1 and modulo p = 2147473409,
 * then an error (SOLVE_ERR_GCD) is reported. This test is not necessary
 * for the computation itself, but fulfilling it implies that G will
 * be recoverable later on from f, g and F. Only a very small proportion
 * of possible polynomials f are not invertible modulo X^n+1 and p.
 *
 * RAM USAGE: 6*n words
 */
#define solve_NTRU   Zn(solve_NTRU)
int solve_NTRU(const ntru_profile *prof, unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g, uint32_t *tmp);

/*
 * Recompute G from f, g and F (using the NTRU equation f*G - g*F = q).
 * This may fail if f is not invertible modulo X^n+1 and modulo
 * p = 2147473409. However, the rest of this key pair generator code takes
 * care never to generate keys with such polynomials f.
 *
 * G is returned as the first n bytes of tmp.
 *
 * Returned value is 1 on success, 0 on error. An error is reported if
 * f is not invertible, or if any of the reconstructed G coefficients is
 * not in the [-lim..+lim] range.
 *
 * RAM USAGE: 3*n words
 */
#define recover_G   Zn(recover_G)
int recover_G(unsigned logn, int32_t q, uint32_t lim,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, uint32_t *restrict tmp);

/* ==================================================================== */

#if NTRUGEN_STATS
/*
 * Statistics are gathered in some global variables. These variables are
 * not thread-safe; thus, this should be disabled by default.
 */
#define stats_hawk_ctt_attempt   Zn(stats_hawk_ctt_attempt)
extern uint32_t stats_hawk_ctt_attempt;
#define stats_hawk_ctt_reject   Zn(stats_hawk_ctt_reject)
extern uint32_t stats_hawk_ctt_reject;
#define stats_solve_attempt   Zn(stats_solve_attempt)
extern uint32_t stats_solve_attempt;
#define stats_solve_err_gcd   Zn(stats_solve_err_gcd)
extern uint32_t stats_solve_err_gcd;
#define stats_solve_err_reduce   Zn(stats_solve_err_reduce)
extern uint32_t stats_solve_err_reduce;
#define stats_solve_err_limit   Zn(stats_solve_err_limit)
extern uint32_t stats_solve_err_limit;
#define stats_solve_success   Zn(stats_solve_success)
extern uint32_t stats_solve_success;
#define stats_compute_w_attempt   Zn(stats_compute_w_attempt)
extern uint32_t stats_compute_w_attempt;
#define stats_compute_w_err_lim1   Zn(stats_compute_w_err_lim1)
extern uint32_t stats_compute_w_err_lim1;
#define stats_compute_w_err_lim2   Zn(stats_compute_w_err_lim2)
extern uint32_t stats_compute_w_err_lim2;
#define stats_compute_w_err_lim3   Zn(stats_compute_w_err_lim3)
extern uint32_t stats_compute_w_err_lim3;
#define stats_compute_w_err_norm   Zn(stats_compute_w_err_norm)
extern uint32_t stats_compute_w_err_norm;
#define stats_compute_w_success   Zn(stats_compute_w_success)
extern uint32_t stats_compute_w_success;

/*
 * Initialize (or reset) the statistics-gathering counters.
 */
#define stats_init   Zn(stats_init)
void stats_init(void);
#endif

/* ==================================================================== */

#endif
