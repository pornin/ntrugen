#include "ng_inner.h"

static const ntru_profile SOLVE_Hawk_256 = {
	1,
	8, 8,
	{ 1, 1, 1, 2, 3, 5, 9, 17, 34, 0, 0 },
	{ 1, 1, 2, 4, 7, 13, 26, 50, 0, 0 },
	{ 1, 1, 1, 2, 3, 3, 3, 4, 0, 0 },
	20,
	{ 0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};

static const ntru_profile SOLVE_Hawk_512 = {
	1,
	9, 9,
	{ 1, 1, 1, 2, 3, 6, 11, 21, 41, 82, 0 },
	{ 1, 2, 3, 5, 8, 16, 31, 61, 121, 0 },
	{ 1, 1, 1, 2, 2, 3, 3, 4, 6, 0 },
	15,
	{ 0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};

static const ntru_profile SOLVE_Hawk_1024 = {
	1,
	10, 10,
	{ 1, 1, 2, 2, 4, 7, 13, 25, 48, 96, 191 },
	{ 1, 2, 3, 5, 10, 19, 37, 72, 143, 284 },
	{ 1, 1, 2, 2, 3, 3, 3, 4, 4, 7 },
	12,
	{ 0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};

#if 0 /* obsolete */
/*
 * Tables for a Gaussian distribution of (f,g), scaled by 32767.
 */

/* Hawk, q = 1, n = 256 -> kmax = 4 */
static const uint16_t gauss_Hawk_256[] = {
	4,
	   16,   305,  2580, 10442, 22325, 30187, 32462, 32751
};

/* Hawk, q = 1, n = 512 -> kmax = 6 */
static const uint16_t gauss_Hawk_512[] = {
	6,
	    3,    37,   286,  1465,  5048, 12026, 20741, 27719,
	31302, 32481, 32730, 32764
};

/* Hawk, q = 1, n = 1024 -> kmax = 8 */
static const uint16_t gauss_Hawk_1024[] = {
	8,
	    2,    17,    89,   377,  1261,  3383,  7347, 13115,
	19652, 25420, 29384, 31506, 32390, 32678, 32750, 32765
};
#endif

TARGET_AVX2
static void
regen_fg_8(int8_t *restrict f, int8_t *restrict g, const void *seed)
{
	size_t seed_len = 16;
#if NTRUGEN_AVX2
	/*
	 * Initialize the four SHAKE256 instances, to run them in parallel.
	 */
	shake_context sc[4];
	for (int i = 0; i < 4; i ++) {
		shake_init(&sc[i], 256);
		shake_inject(&sc[i], seed, seed_len);
		uint8_t ix = (uint8_t)i;
		shake_inject(&sc[i], &ix, 1);
	}
	shake_x4_context scx4;
	shake_x4_flip(&scx4, sc);
	__m256i ym4 = _mm256_set1_epi8(0x0F);
	__m256i ytt = _mm256_setr_epi8(
		-2, -1, -1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 1, 1, 2,
		-2, -1, -1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 1, 1, 2);
	for (size_t u = 0; u < 512; u += 64) {
		union {
			__m256i y;
			uint64_t q[4];
		} buf;
		shake_x4_extract_words(&scx4, buf.q, 1);
		__m256i y0 = _mm256_and_si256(buf.y, ym4);
		__m256i y1 = _mm256_and_si256(_mm256_srli_epi16(buf.y, 4), ym4);
		y0 = _mm256_shuffle_epi8(ytt, y0);
		y1 = _mm256_shuffle_epi8(ytt, y1);
		__m256i yv0 = _mm256_unpacklo_epi8(y0, y1);
		__m256i yv1 = _mm256_unpackhi_epi8(y0, y1);
		__m256i yd0 = _mm256_permute2x128_si256(yv0, yv1, 0x20);
		__m256i yd1 = _mm256_permute2x128_si256(yv0, yv1, 0x31);
		if (u < 256) {
			_mm256_storeu_si256((__m256i *)(f + u +   0), yd0);
			_mm256_storeu_si256((__m256i *)(f + u +  32), yd1);
		} else {
			_mm256_storeu_si256((__m256i *)(g + u - 256), yd0);
			_mm256_storeu_si256((__m256i *)(g + u - 224), yd1);
		}
	}
#else // NTRUGEN_AVX2
	for (size_t j = 0; j < 4; j ++) {
		shake_context sc;
		shake_init(&sc, 256);
		shake_inject(&sc, seed, seed_len);
		uint8_t jx = (uint8_t)j;
		shake_inject(&sc, &jx, 1);
		shake_flip(&sc);
		for (size_t u = 0; u < 512; u += 64) {
			uint8_t qb[8];
			shake_extract(&sc, qb, 8);
			uint64_t q = dec64le(qb);
			q = (q & (uint64_t)0x5555555555555555)
				+ ((q >> 1) & (uint64_t)0x5555555555555555);
			q = (q & (uint64_t)0x3333333333333333)
				+ ((q >> 2) & (uint64_t)0x3333333333333333);
			int8_t vv[16];
			for (int i = 0; i < 16; i ++) {
				vv[i] = (int)(q & 0x0F) - 2;
				q >>= 4;
			}
			if (u < 256) {
				memcpy(f + u + (j << 4), vv, 16);
			} else {
				memcpy(g + (u - 256) + (j << 4), vv, 16);
			}
		}
	}
#endif // NTRUGEN_AVX2
}

TARGET_AVX2
static void
regen_fg_9(int8_t *restrict f, int8_t *restrict g, const void *seed)
{
	size_t seed_len = 24;
#if NTRUGEN_AVX2
	/*
	 * Initialize the four SHAKE256 instances, to run them in parallel.
	 */
	shake_context sc[4];
	for (int i = 0; i < 4; i ++) {
		shake_init(&sc[i], 256);
		shake_inject(&sc[i], seed, seed_len);
		uint8_t ix = (uint8_t)i;
		shake_inject(&sc[i], &ix, 1);
	}
	shake_x4_context scx4;
	shake_x4_flip(&scx4, sc);
	__m256i ym4 = _mm256_set1_epi8(0x0F);
	__m256i ytt = _mm256_setr_epi8(
		-2, -1, -1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 1, 1, 2,
		-2, -1, -1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 1, 1, 2);
	for (size_t u = 0; u < 1024; u += 32) {
		union {
			__m256i y;
			uint64_t q[4];
		} buf;
		shake_x4_extract_words(&scx4, buf.q, 1);
		__m256i y0 = _mm256_and_si256(buf.y, ym4);
		__m256i y1 = _mm256_and_si256(_mm256_srli_epi16(buf.y, 4), ym4);
		y0 = _mm256_shuffle_epi8(ytt, y0);
		y1 = _mm256_shuffle_epi8(ytt, y1);
		__m256i yv = _mm256_add_epi8(y0, y1);
		if (u < 512) {
			_mm256_storeu_si256((__m256i *)(f + u +   0), yv);
		} else {
			_mm256_storeu_si256((__m256i *)(g + u - 512), yv);
		}
	}
#else // NTRUGEN_AVX2
	for (size_t j = 0; j < 4; j ++) {
		shake_context sc;
		shake_init(&sc, 256);
		shake_inject(&sc, seed, seed_len);
		uint8_t jx = (uint8_t)j;
		shake_inject(&sc, &jx, 1);
		shake_flip(&sc);
		for (size_t u = 0; u < 1024; u += 32) {
			uint8_t qb[8];
			shake_extract(&sc, qb, 8);
			uint64_t q = dec64le(qb);
			q = (q & (uint64_t)0x5555555555555555)
				+ ((q >> 1) & (uint64_t)0x5555555555555555);
			q = (q & (uint64_t)0x3333333333333333)
				+ ((q >> 2) & (uint64_t)0x3333333333333333);
			q = (q & (uint64_t)0x0F0F0F0F0F0F0F0F)
				+ ((q >> 4) & (uint64_t)0x0F0F0F0F0F0F0F0F);
			int8_t vv[8];
			for (int i = 0; i < 8; i ++) {
				vv[i] = (int)(q & 0xFF) - 4;
				q >>= 8;
			}
			if (u < 512) {
				memcpy(f + u + (j << 3), vv, 8);
			} else {
				memcpy(g + (u - 512) + (j << 3), vv, 8);
			}
		}
	}
#endif // NTRUGEN_AVX2
}

TARGET_AVX2
static void
regen_fg_10(int8_t *restrict f, int8_t *restrict g, const void *seed)
{
	size_t seed_len = 40;
#if NTRUGEN_AVX2
	/*
	 * Initialize the four SHAKE256 instances, to run them in parallel.
	 */
	shake_context sc[4];
	for (int i = 0; i < 4; i ++) {
		shake_init(&sc[i], 256);
		shake_inject(&sc[i], seed, seed_len);
		uint8_t ix = (uint8_t)i;
		shake_inject(&sc[i], &ix, 1);
	}
	shake_x4_context scx4;
	shake_x4_flip(&scx4, sc);
	__m256i ym4 = _mm256_set1_epi8(0x0F);
	__m256i ytt = _mm256_setr_epi8(
		-2, -1, -1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 1, 1, 2,
		-2, -1, -1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 1, 1, 2);
	__m256i ys = _mm256_setr_epi8(
		0, 2, 4, 6, 8, 10, 12, 14,
		-1, -1, -1, -1, -1, -1, -1, -1,
		0, 2, 4, 6, 8, 10, 12, 14,
		-1, -1, -1, -1, -1, -1, -1, -1);
	for (size_t u = 0; u < 2048; u += 16) {
		union {
			__m256i y;
			uint64_t q[4];
		} buf;
		shake_x4_extract_words(&scx4, buf.q, 1);
		__m256i y0 = _mm256_and_si256(buf.y, ym4);
		__m256i y1 = _mm256_and_si256(_mm256_srli_epi16(buf.y, 4), ym4);
		y0 = _mm256_shuffle_epi8(ytt, y0);
		y1 = _mm256_shuffle_epi8(ytt, y1);
		__m256i yv = _mm256_add_epi8(y0, y1);
		yv = _mm256_add_epi8(yv, _mm256_srli_epi16(yv, 8));
		yv = _mm256_shuffle_epi8(yv, ys);
		yv = _mm256_permute4x64_epi64(yv, 0xD8);
		__m128i xv = _mm256_castsi256_si128(yv);
		if (u < 1024) {
			_mm_storeu_si128((__m128i *)(f + u +    0), xv);
		} else {
			_mm_storeu_si128((__m128i *)(g + u - 1024), xv);
		}
	}
#else // NTRUGEN_AVX2
	for (size_t j = 0; j < 4; j ++) {
		shake_context sc;
		shake_init(&sc, 256);
		shake_inject(&sc, seed, seed_len);
		uint8_t jx = (uint8_t)j;
		shake_inject(&sc, &jx, 1);
		shake_flip(&sc);
		for (size_t u = 0; u < 2048; u += 16) {
			uint8_t qb[8];
			shake_extract(&sc, qb, 8);
			uint64_t q = dec64le(qb);
			q = (q & (uint64_t)0x5555555555555555)
				+ ((q >> 1) & (uint64_t)0x5555555555555555);
			q = (q & (uint64_t)0x3333333333333333)
				+ ((q >> 2) & (uint64_t)0x3333333333333333);
			q = (q & (uint64_t)0x0F0F0F0F0F0F0F0F)
				+ ((q >> 4) & (uint64_t)0x0F0F0F0F0F0F0F0F);
			q = (q & (uint64_t)0x00FF00FF00FF00FF)
				+ ((q >> 8) & (uint64_t)0x00FF00FF00FF00FF);
			int8_t vv[4];
			for (int i = 0; i < 4; i ++) {
				vv[i] = (int)(q & 0xFFFF) - 8;
				q >>= 16;
			}
			if (u < 1024) {
				memcpy(f + u + (j << 2), vv, 4);
			} else {
				memcpy(g + (u - 1024) + (j << 2), vv, 4);
			}
		}
	}
#endif // NTRUGEN_AVX2
}

/* see ntrugen.h */
TARGET_AVX2
void
Hawk_regen_fg(unsigned logn,
	int8_t *restrict f, int8_t *restrict g, const void *seed)
{
	switch (logn) {
	case 8:
		regen_fg_8(f, g, seed);
		break;
	case 9:
		regen_fg_9(f, g, seed);
		break;
	default:
		regen_fg_10(f, g, seed);
		break;
	}

#if 0
	const uint16_t *tab;
	switch (logn) {
	case 8:   tab = gauss_Hawk_256;   break;
	case 9:   tab = gauss_Hawk_512;   break;
	default:  tab = gauss_Hawk_1024;  break;
	}
	size_t n = (size_t)1 << logn;
	size_t kmax = tab[0];
	size_t seed_len = 8 + ((size_t)1 << (logn - 5));

#if NTRUGEN_AVX2

	/*
	 * Initialize the four SHAKE256 instances, to run them in parallel.
	 */
	shake_context sc[4];
	for (int i = 0; i < 4; i ++) {
		shake_init(&sc[i], 256);
		shake_inject(&sc[i], seed, seed_len);
		uint8_t ix = (uint8_t)i;
		shake_inject(&sc[i], &ix, 1);
	}
	shake_x4_context scx4;
	shake_x4_flip(&scx4, sc);

	__m256i ytab[16];
	for (size_t u = 0; u < (kmax << 1); u ++) {
		ytab[u] = _mm256_set1_epi16(tab[u + 1]);
	}
	__m256i yb = _mm256_set1_epi16(-(int)kmax);
	__m256i ys = _mm256_setr_epi8(
		0, 2, 4, 6, 8, 10, 12, 14,
		-1, -1, -1, -1, -1, -1, -1, -1,
		0, 2, 4, 6, 8, 10, 12, 14,
		-1, -1, -1, -1, -1, -1, -1, -1);
	__m256i y15 = _mm256_set1_epi16(0x7FFF);

	for (size_t u = 0; u < (n << 1); u += 16) {
		union {
			__m256i y;
			uint64_t q[4];
		} buf;
		shake_x4_extract_words(&scx4, buf.q, 1);
		__m256i yr = _mm256_and_si256(buf.y, y15);
		__m256i yv = yb;
		for (size_t k = 0; k < (kmax << 1); k ++) {
			yv = _mm256_sub_epi16(yv,
				_mm256_cmpgt_epi16(yr, ytab[k]));
		}
		yv = _mm256_shuffle_epi8(yv, ys);
		yv = _mm256_permute4x64_epi64(yv, 0xD8);
		__m128i xv = _mm256_castsi256_si128(yv);
		if (u < n) {
			_mm_storeu_si128((__m128i *)(f + u), xv);
		} else {
			_mm_storeu_si128((__m128i *)(g + (u - n)), xv);
		}
	}

#else // NTRUGEN_AVX2

	for (size_t j = 0; j < 4; j ++) {
		shake_context sc;
		shake_init(&sc, 256);
		shake_inject(&sc, seed, seed_len);
		uint8_t jx = (uint8_t)j;
		shake_inject(&sc, &jx, 1);
		shake_flip(&sc);
		for (size_t u = 0; u < (n << 1); u += 16) {
			uint8_t qb[8];
			shake_extract(&sc, qb, 8);
			uint64_t q = dec64le(qb);
			int8_t vv[4];
			for (int i = 0; i < 4; i ++) {
				uint32_t x = (uint32_t)q & 0x7FFF;
				q >>= 16;
				uint32_t v = -(uint32_t)kmax;
				for (size_t k = 1; k <= (kmax << 1); k ++) {
					v += ((uint32_t)tab[k] - x) >> 31;
				}
				vv[i] = (int8_t)*(int32_t *)&v;
			}
			if (u < n) {
				memcpy(f + u + (j << 2), vv, 4);
			} else {
				memcpy(g + (u - n) + (j << 2), vv, 4);
			}
		}
	}

#endif // NTRUGEN_AVX2
#endif
}

TARGET_AVX2
static unsigned
parity(unsigned logn, int8_t *f)
{
	size_t n = (size_t)1 << logn;
#if NTRUGEN_AVX2
	__m256i yr = _mm256_setzero_si256();
	for (size_t u = 0; u < n; u += 32) {
		__m256i y = _mm256_loadu_si256((const __m256i *)(f + u));
		yr = _mm256_xor_si256(y, yr);
	}
	uint32_t r = (uint32_t)_mm256_movemask_epi8(_mm256_slli_epi16(yr, 7));
	r ^= (r >> 16);
	r ^= (r >> 8);
	r ^= (r >> 4);
	r ^= (r >> 2);
	r ^= (r >> 1);
	return (unsigned)(r & 1);
#else // NTRUGEN_AVX2
	unsigned pp = 0;
	for (size_t u = 0; u < n; u ++) {
		pp += *(uint8_t *)&f[u];
	}
	return pp & 1;
#endif // NTRUGEN_AVX2
}

/*
 * Limits for q00, q01 and q11 (maximum bit size of the absolute value of
 * a coefficient, excluding q00[0] and q11[0]).
 */
static const int8_t bits_lim00[11] = {
	0, 0, 0, 0, 0, 0, 0, 0,  9,  9, 10
};
static const int8_t bits_lim01[11] = {
	0, 0, 0, 0, 0, 0, 0, 0, 11, 12, 14
};
static const int8_t bits_lim11[11] = {
	0, 0, 0, 0, 0, 0, 0, 0, 13, 15, 17
};

/*
 * Given f, g, F and G, compute:
 *   q00 = f*adj(f) + g*adj(g)
 *   q01 = F*adj(f) + G*adj(g)
 *   q11 = F*adj(f) + G*adj(G)
 * The three polynomials are stored at the start of tmp[], in that order:
 *    q00   int16_t[n]
 *    q01   int16_t[n]
 *    q11   int32_t[n]
 * Note: q00 and q11 are auto-adjoint.
 *
 * Return value is 1 on success, 0 on error. An error is reported if
 * any of the coefficients does not comply with the following limits:
 *    -32768 <= q00[0] < 32768
 *    -lim00 <= q00[u] < +lim00    for u = 1 to n  
 *    -lim01 <= q01[u] < +lim01    for u = 0 to n  
 *    -lim11 <= q11[u] < +lim11    for u = 1 to n
 * An error is also reported if q00 turns out not to be invertible modulo
 * X^n+1 and modulo 2147473409.
 *
 * RAM USAGE: 5*n words
 */
static int
make_q001(unsigned logn,
	int lim00, int lim01, int32_t lim11,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	uint32_t *restrict tmp)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	uint32_t R2 = PRIMES[0].R2;

	uint32_t *t1 = tmp;
	uint32_t *t2 = t1 + n;
	uint32_t *t3 = t2 + n;
	uint32_t *t4 = t3 + n;
	uint32_t *t5 = t4 + n;

	mp_mkgm(logn, t1, PRIMES[0].g, p, p0i);
	poly_mp_set_small(logn, t2, f, p);
	poly_mp_set_small(logn, t3, g, p);
	poly_mp_set_small(logn, t4, F, p);
	poly_mp_set_small(logn, t5, G, p);
	mp_NTT(logn, t2, t1, p, p0i);
	mp_NTT(logn, t3, t1, p, p0i);
	mp_NTT(logn, t4, t1, p, p0i);
	mp_NTT(logn, t5, t1, p, p0i);
	for (size_t u = 0; u < hn; u ++) {
		uint32_t xf = t2[u];
		uint32_t xfa = t2[(n - 1) - u];
		uint32_t xg = t3[u];
		uint32_t xga = t3[(n - 1) - u];
		uint32_t xF = t4[u];
		uint32_t xFa = t4[(n - 1) - u];
		uint32_t xG = t5[u];
		uint32_t xGa = t5[(n - 1) - u];

		uint32_t xq00 = mp_montymul(R2, mp_add(
			mp_montymul(xf, xfa, p, p0i),
			mp_montymul(xg, xga, p, p0i), p), p, p0i);
		uint32_t xq11 = mp_montymul(R2, mp_add(
			mp_montymul(xF, xFa, p, p0i),
			mp_montymul(xG, xGa, p, p0i), p), p, p0i);
		uint32_t xq01_0 = mp_montymul(R2, mp_add(
			mp_montymul(xF, xfa, p, p0i),
			mp_montymul(xG, xga, p, p0i), p), p, p0i);
		uint32_t xq01_1 = mp_montymul(R2, mp_add(
			mp_montymul(xFa, xf, p, p0i),
			mp_montymul(xGa, xg, p, p0i), p), p, p0i);
		if (xq00 == 0) {
			/* q00 is not invertible mod X^n+1 mod p,
			   we report the error right away (key is invalid
			   and will be discarded). */
			return 0;
		}
		t3[u] = xq00;
		t3[(n - 1) - u] = xq00;
		t4[u] = xq01_0;
		t4[(n - 1) - u] = xq01_1;
		t5[u] = xq11;
		t5[(n - 1) - u] = xq11;
	}
	mp_mkigm(logn, t1, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t3, t1, p, p0i);
	mp_iNTT(logn, t4, t1, p, p0i);
	mp_iNTT(logn, t5, t1, p, p0i);

	/*
	 * t1 and t2 are free, we reuse them for the output.
	 */
	int16_t *q00 = (int16_t *)tmp;
	int16_t *q01 = q00 + n;
	int32_t *q11 = (int32_t *)(q01 + n);
	for (size_t u = 0; u < n; u ++) {
		int32_t xq00 = mp_norm(t3[u], p);
		int32_t xq01 = mp_norm(t4[u], p);
		int32_t xq11 = mp_norm(t5[u], p);
		if (u == 0) {
			if (xq00 < -32768 || xq00 > +32767) {
				return 0;
			}
		} else {
			if (xq00 <= -lim00 || xq00 >= +lim00) {
				return 0;
			}
			if (xq11 <= -lim11 || xq11 >= +lim11) {
				return 0;
			}
		}
		if (xq01 <= -lim01 || xq01 >= +lim01) {
			return 0;
		}
		q00[u] = (int16_t)xq00;
		q01[u] = (int16_t)xq01;
		q11[u] = xq11;
	}

	return 1;
}

/* see ntrugen.h */
int
Hawk_keygen(unsigned logn,
	int8_t *restrict f, int8_t *restrict g,
	int8_t *restrict F, int8_t *restrict G,
	int16_t *restrict q00, int16_t *restrict q01, int32_t *restrict q11,
	void *seed, ntrugen_rng rng, void *restrict rng_context,
	void *restrict tmp, size_t tmp_len)
{
	/*
	 * Ensure that the tmp[] buffer has proper alignment for 64-bit
	 * access.
	 */
	if (tmp_len < 7) {
		return -1;
	}
	if (logn < 2 || logn > 10) {
		return -1;
	}
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 7) & ~(uintptr_t)7;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint32_t *tt32 = (void *)utmp2;
	if (tmp_len < ((size_t)24 << logn)) {
		return -1;
	}

	uint32_t l2low;
	fxr d0high;
	const ntru_profile *prof;
	switch (logn) {
	case 8:
		l2low = 556;
		d0high = fxr_of_scaled32(17179869); /* 1/250 */
		prof = &SOLVE_Hawk_256;
		break;
	case 9:
		l2low = 2080;
		d0high = fxr_of_scaled32(4294967); /* 1/1000 */
		prof = &SOLVE_Hawk_512;
		break;
	case 10:
		l2low = 7981;
		d0high = fxr_of_scaled32(1431655); /* 1/3000 */
		prof = &SOLVE_Hawk_1024;
		break;
	default:
		/*
		 * Other degrees are not supported.
		 */
		return -1;
	}
	int lim00 = 1 << bits_lim00[logn];
	int lim01 = 1 << bits_lim01[logn];
	int lim11 = (int32_t)1 << bits_lim11[logn];

	uint8_t seed_buf[40];
	size_t seed_len = 8 + ((size_t)1 << (logn - 5));

	for (;;) {
		/*
		 * Generate f and g.
		 */
		rng(rng_context, seed_buf, seed_len);
		Hawk_regen_fg(logn, f, g, seed_buf);

		/*
		 * Start again if f and g are not both odd.
		 */
		if (parity(logn, f) != 1 || parity(logn, g) != 1) {
			continue;
		}

		/*
		 * Check that (f,g) has an acceptable norm; this is a
		 * _minimum_ bound (2*n*sigma_sec^2).
		 */
		uint32_t norm2_fg = poly_sqnorm(logn, f) + poly_sqnorm(logn, g);
		if (norm2_fg < l2low) {
			continue;
		}

		/*
		 * Check that f*adj(f) + g*adj(g) is invertible modulo
		 * X^n+1 mod p1 (with p1 = 2147473409 = PRIMES[0].p).
		 * We also output f*adj(f) + g*adj(g) into t1.
		 */
		int invertible = 1;
		size_t n = (size_t)1 << logn;
		size_t hn = n >> 1;
		uint32_t *t1 = tt32;
		uint32_t *t2 = t1 + n;
		uint32_t *t3 = t2 + n;
		uint32_t *t4 = t3 + n;
		uint32_t p = PRIMES[0].p;
		uint32_t p0i = PRIMES[0].p0i;
		uint32_t R2 = PRIMES[0].R2;
		mp_mkgmigm(logn, t1, t2, PRIMES[0].g, PRIMES[0].ig, p, p0i);
		for (size_t u = 0; u < n; u ++) {
			t3[u] = mp_set(f[u], p);
			t4[u] = mp_set(g[u], p);
		}
		mp_NTT(logn, t3, t1, p, p0i);
		mp_NTT(logn, t4, t1, p, p0i);
		for (size_t u = 0; u < n; u ++) {
			uint32_t x = mp_add(
				mp_montymul(t3[u], t3[(n - 1) - u], p, p0i),
				mp_montymul(t4[u], t4[(n - 1) - u], p, p0i), p);
			/*
			 * Value x is in anti-Montgomery representation;
			 * it is enough to test invertibility, but we need
			 * the actual value for the next test on (f,g).
			 */
			if (x == 0) {
				invertible = 0;
				break;
			}
			x = mp_montymul(R2, x, p, p0i);
			t1[u] = x;
		}
		if (!invertible) {
			continue;
		}
		/* Get the plain f*adj(f) + g*adj(g) */
		mp_iNTT(logn, t1, t2, p, p0i);
		for (size_t u = 0; u < n; u ++) {
			t1[u] = (uint32_t)mp_norm(t1[u], p);
		}

		/*
		 * Also check that f*adj(f) + g*adj(g) is invertible modulo
		 * X^n+1 mod p2 (with p2 = 2147389441 = PRIMES[1].p).
		 */
		p = PRIMES[1].p;
		p0i = PRIMES[1].p0i;
		for (size_t u = 0; u < n; u ++) {
			t2[u] = mp_set(*(int32_t *)&t1[u], p);
		}
		mp_mkgm(logn, t3, PRIMES[1].g, p, p0i);
		mp_NTT(logn, t2, t3, p, p0i);
		for (size_t u = 0; u < n; u ++) {
			if (t2[u] == 0) {
				invertible = 0;
				break;
			}
		}
		if (!invertible) {
			continue;
		}

		/*
		 * Check that the constant term of 1/(f*adj(f) + g*adj(g))
		 * is small enough.
		 */
#if NTRUGEN_STATS
		stats_hawk_ctt_attempt ++;
#endif
		fxr *rt1 = (fxr *)t2;
		for (size_t u = 0; u < n; u ++) {
			rt1[u] = fxr_of(*(int32_t *)&t1[u]);
		}
		vect_FFT(logn, rt1);
		for (size_t u = 0; u < hn; u ++) {
			rt1[u] = fxr_inv(rt1[u]);
		}
		/* Normally the values are already zero, or close to zero
		   in case of loss of precision. We force them to zero. */
		for (size_t u = hn; u < n; u ++) {
			rt1[u] = fxr_zero;
		}
		vect_iFFT(logn, rt1);

		if (fxr_lt(d0high, rt1[0])) {
#if NTRUGEN_STATS
			stats_hawk_ctt_reject ++;
#endif
			continue;
		}

		/*
		 * Solve the NTRU equation.
		 */
#if NTRUGEN_STATS
		stats_solve_attempt ++;
#endif
		int err = solve_NTRU(prof, logn, f, g, tt32);
		switch (err) {
		case SOLVE_OK:
#if NTRUGEN_STATS
			stats_solve_success ++;
#endif
			break;
#if NTRUGEN_STATS
		case SOLVE_ERR_GCD:
			stats_solve_err_gcd ++;
			continue;
		case SOLVE_ERR_REDUCE:
			stats_solve_err_reduce ++;
			continue;
		case SOLVE_ERR_LIMIT:
			stats_solve_err_limit ++;
			continue;
#endif
		default:
			continue;
		}

		/*
		 * F and G are at the start of tt32[].
		 */
		int8_t *tF = (int8_t *)tt32;
		int8_t *tG = tF + n;

		/*
		 * Compute q00, q01 and q1, and check that they are in the
		 * expected range.
		 *
		 * F and G use the first 2*n bytes = hn words.
		 */
		if (!make_q001(logn, lim00, lim01, lim11,
			f, g, tF, tG, (uint32_t *)(tG + n)))
		{
#if NTRUGEN_STATS
			stats_solve_err_limit ++;
#endif
			continue;
		}

		int16_t *tq00 = (int16_t *)(tG + n);
		int16_t *tq01 = tq00 + n;
		int32_t *tq11 = (int32_t *)(tq01 + n);
		uint8_t *tseed = (uint8_t *)(tq11 + n);
		memmove(tseed, seed_buf, seed_len);

		/*
		 * Return the computed F, G, q00, q01, q11 and seed.
		 */
		if (F != NULL) {
			memmove(F, tF, n);
		}
		if (G != NULL) {
			memmove(G, tG, n);
		}
		if (q00 != NULL) {
			memmove(q00, tq00, n * sizeof *tq00);
		}
		if (q01 != NULL) {
			memmove(q01, tq01, n * sizeof *tq01);
		}
		if (q11 != NULL) {
			memmove(q11, tq11, n * sizeof *tq11);
		}
		if (seed != NULL) {
			memmove(seed, tseed, seed_len);
		}
		if (tt32 != tmp) {
			memmove(tmp, tt32, 10 * n + seed_len);
		}

		return 0;
	}
}

/* see ntrugen.h */
int
Hawk_recover_G(unsigned logn,
	int8_t *restrict G,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, void *restrict tmp, size_t tmp_len)
{
	const ntru_profile *prof;
	switch (logn) {
	case 8:   prof = &SOLVE_Hawk_256;   break;
	case 9:   prof = &SOLVE_Hawk_512;   break;
	case 10:  prof = &SOLVE_Hawk_1024;  break;
	default:
		return -1;
	}

	/*
	 * Ensure that the tmp[] buffer has proper alignment for 64-bit
	 * access.
	 */
	if (tmp_len < 7) {
		return -1;
	}
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 7) & ~(uintptr_t)7;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint32_t *tt32 = (void *)utmp2;
	if (tmp_len < ((size_t)12 << logn)) {
		return -1;
	}

	int r = recover_G(logn,
		(int32_t)prof->q, prof->coeff_FG_limit[logn],
		f, g, F, tt32);
	size_t n = (size_t)1 << logn;
	if (G != NULL) {
		memmove(G, tt32, n);
	}
	if ((void *)tt32 != tmp) {
		memmove(tmp, tt32, n);
	}

	return r - 1;
}

/* see ntrugen.h */
int
Hawk_recover_qq(unsigned logn,
	int16_t *restrict q00, int16_t *restrict q01, int32_t *restrict q11,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	void *restrict tmp, size_t tmp_len)
{
	int lim00, lim01;
	int32_t lim11;
	switch (logn) {
	case 8:
	case 9:
	case 10:
		lim00 = 1 << bits_lim00[logn];
		lim01 = 1 << bits_lim01[logn];
		lim11 = (int32_t)1 << bits_lim11[logn];
		break;
	default:
		return -1;
	}

	/*
	 * Ensure that the tmp[] buffer has proper alignment for 64-bit
	 * access.
	 */
	if (tmp_len < 7) {
		return -1;
	}
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 7) & ~(uintptr_t)7;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint32_t *tt32 = (void *)utmp2;
	if (tmp_len < ((size_t)20 << logn)) {
		return -1;
	}

	if (!make_q001(logn, lim00, lim01, lim11, f, g, F, G, tt32)) {
		return -1;
	}
	size_t n = (size_t)1 << logn;
	int16_t *tq00 = (int16_t *)tt32;
	int16_t *tq01 = tq00 + n;
	int32_t *tq11 = (int32_t *)(tq01 + n);
	if (q00 != NULL) {
		memmove(q00, tq00, n * sizeof *tq00);
	}
	if (q01 != NULL) {
		memmove(q01, tq01, n * sizeof *tq01);
	}
	if (q11 != NULL) {
		memmove(q11, tq11, n * sizeof *tq11);
	}
	if ((void *)tt32 != tmp) {
		memmove(tmp, tt32, n * 8);
	}
	return 0;
}
