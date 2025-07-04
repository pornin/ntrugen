#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "ntrugen.h"
#include "ntrugen_prng.h"

/*
 * We test some inner functions, so we use the inner API.
 */
#include "ng_inner.h"

static void *
xmalloc(size_t len)
{
	void *buf;

	if (len == 0) {
		return NULL;
	}
	buf = malloc(len);
	if (buf == NULL) {
		fprintf(stderr, "memory allocation error\n");
		exit(EXIT_FAILURE);
	}
	return buf;
}

static void
xfree(void *buf)
{
	if (buf != NULL) {
		free(buf);
	}
}

static size_t
hextobin(uint8_t *buf, size_t max_len, const char *src)
{
	size_t u;
	int acc, z;

	u = 0;
	acc = 0;
	z = 0;
	for (;;) {
		int c;

		c = *src ++;
		if (c == 0) {
			if (z) {
				fprintf(stderr, "Lone hex nibble\n");
				exit(EXIT_FAILURE);
			}
			return u;
		}
		if (c >= '0' && c <= '9') {
			c -= '0';
		} else if (c >= 'A' && c <= 'F') {
			c -= 'A' - 10;
		} else if (c >= 'a' && c <= 'f') {
			c -= 'a' - 10;
		} else if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
			continue;
		} else {
			fprintf(stderr, "Not an hex digit: U+%04X\n",
				(unsigned)c);
			exit(EXIT_FAILURE);
		}
		if (z) {
			if (u >= max_len) {
				fprintf(stderr,
					"Hex string too long for buffer\n");
				exit(EXIT_FAILURE);
			}
			buf[u ++] = (unsigned char)((acc << 4) + c);
		} else {
			acc = c;
		}
		z = !z;
	}
}

#define HEXTOBIN(dst, src)   do { \
		if (hextobin(dst, sizeof(dst), src) != sizeof(dst)) { \
			fprintf(stderr, "Wrong hexdec length\n"); \
			exit(EXIT_FAILURE); \
		} \
	} while (0)

static void
check_eq(const void *a, const void *b, size_t len, const char *banner)
{
	size_t u;

	if (memcmp(a, b, len) == 0) {
		return;
	}
	fprintf(stderr, "%s: wrong value:\n", banner);
	fprintf(stderr, "a: ");
	for (u = 0; u < len; u ++) {
		fprintf(stderr, "%02x", ((const unsigned char *)a)[u]);
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "b: ");
	for (u = 0; u < len; u ++) {
		fprintf(stderr, "%02x", ((const unsigned char *)b)[u]);
	}
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}

/*
 * A perfunctory SHA-1 implementation, for aggregate reproducible tests.
 */
static void
sha1_round_inner(const uint8_t *buf, uint32_t *val)
{
#define F(B, C, D)     ((((C) ^ (D)) & (B)) ^ (D))
#define G(B, C, D)     ((B) ^ (C) ^ (D))
#define H(B, C, D)     (((D) & (C)) | (((D) | (C)) & (B)))
#define I(B, C, D)     G(B, C, D)

#define ROTL(x, n)    (((x) << (n)) | ((x) >> (32 - (n))))

#define K1     ((uint32_t)0x5A827999)
#define K2     ((uint32_t)0x6ED9EBA1)
#define K3     ((uint32_t)0x8F1BBCDC)
#define K4     ((uint32_t)0xCA62C1D6)

	uint32_t m[80];
	for (int i = 0; i < 16; i ++) {
		m[i] = ((uint32_t)buf[(i << 2) + 0] << 24)
			| ((uint32_t)buf[(i << 2) + 1] << 16)
			| ((uint32_t)buf[(i << 2) + 2] << 8)
			| (uint32_t)buf[(i << 2) + 3];
	}
	for (int i = 16; i < 80; i ++) {
		uint32_t x = m[i - 3] ^ m[i - 8] ^ m[i - 14] ^ m[i - 16];
		m[i] = ROTL(x, 1);
	}

	uint32_t a = val[0];
	uint32_t b = val[1];
	uint32_t c = val[2];
	uint32_t d = val[3];
	uint32_t e = val[4];
	for (int i = 0; i < 20; i += 5) {
		e += ROTL(a, 5) + F(b, c, d) + K1 + m[i + 0]; b = ROTL(b, 30);
		d += ROTL(e, 5) + F(a, b, c) + K1 + m[i + 1]; a = ROTL(a, 30);
		c += ROTL(d, 5) + F(e, a, b) + K1 + m[i + 2]; e = ROTL(e, 30);
		b += ROTL(c, 5) + F(d, e, a) + K1 + m[i + 3]; d = ROTL(d, 30);
		a += ROTL(b, 5) + F(c, d, e) + K1 + m[i + 4]; c = ROTL(c, 30);
	}
	for (int i = 20; i < 40; i += 5) {
		e += ROTL(a, 5) + G(b, c, d) + K2 + m[i + 0]; b = ROTL(b, 30);
		d += ROTL(e, 5) + G(a, b, c) + K2 + m[i + 1]; a = ROTL(a, 30);
		c += ROTL(d, 5) + G(e, a, b) + K2 + m[i + 2]; e = ROTL(e, 30);
		b += ROTL(c, 5) + G(d, e, a) + K2 + m[i + 3]; d = ROTL(d, 30);
		a += ROTL(b, 5) + G(c, d, e) + K2 + m[i + 4]; c = ROTL(c, 30);
	}
	for (int i = 40; i < 60; i += 5) {
		e += ROTL(a, 5) + H(b, c, d) + K3 + m[i + 0]; b = ROTL(b, 30);
		d += ROTL(e, 5) + H(a, b, c) + K3 + m[i + 1]; a = ROTL(a, 30);
		c += ROTL(d, 5) + H(e, a, b) + K3 + m[i + 2]; e = ROTL(e, 30);
		b += ROTL(c, 5) + H(d, e, a) + K3 + m[i + 3]; d = ROTL(d, 30);
		a += ROTL(b, 5) + H(c, d, e) + K3 + m[i + 4]; c = ROTL(c, 30);
	}
	for (int i = 60; i < 80; i += 5) {
		e += ROTL(a, 5) + I(b, c, d) + K4 + m[i + 0]; b = ROTL(b, 30);
		d += ROTL(e, 5) + I(a, b, c) + K4 + m[i + 1]; a = ROTL(a, 30);
		c += ROTL(d, 5) + I(e, a, b) + K4 + m[i + 2]; e = ROTL(e, 30);
		b += ROTL(c, 5) + I(d, e, a) + K4 + m[i + 3]; d = ROTL(d, 30);
		a += ROTL(b, 5) + I(c, d, e) + K4 + m[i + 4]; c = ROTL(c, 30);
	}

	val[0] += a;
	val[1] += b;
	val[2] += c;
	val[3] += d;
	val[4] += e;

#undef F
#undef G
#undef H
#undef I
#undef ROTL
#undef K1
#undef K2
#undef K3
#undef K4
}

typedef struct {
	uint8_t buf[64];
	uint32_t val[5];
	uint64_t count;
} sha1_context;

static void
sha1_init(sha1_context *sc)
{
	static const uint32_t IV[5] = {
		0x67452301, 0xEFCDAB89, 0x98BADCFE, 0x10325476, 0xC3D2E1F0
	};

	memset(sc->buf, 0, sizeof sc->buf);
	memcpy(sc->val, IV, sizeof sc->val);
	sc->count = 0;
}

static void
sha1_update(sha1_context *sc, const void *data, size_t len)
{
	const uint8_t *buf = data;
	size_t ptr = (size_t)sc->count & 63;
	sc->count += (uint64_t)len;
	while (len > 0) {
		size_t clen = 64 - ptr;
		if (clen > len) {
			clen = len;
		}
		memcpy(sc->buf + ptr, buf, clen);
		buf += clen;
		len -= clen;
		ptr += clen;
		if (ptr == 64) {
			sha1_round_inner(sc->buf, sc->val);
			ptr = 0;
		}
	}
}

static void
sha1_update_u16le(sha1_context *sc, const uint16_t *data, size_t num)
{
	uint8_t tmp[2];
	for (size_t u = 0; u < num; u ++) {
		tmp[0] = (uint8_t)data[u];
		tmp[1] = (uint8_t)(data[u] >> 8);
		sha1_update(sc, tmp, 2);
	}
}

static void
sha1_update_u32le(sha1_context *sc, const uint32_t *data, size_t num)
{
	uint8_t tmp[4];
	for (size_t u = 0; u < num; u ++) {
		tmp[0] = (uint8_t)data[u];
		tmp[1] = (uint8_t)(data[u] >> 8);
		tmp[2] = (uint8_t)(data[u] >> 16);
		tmp[3] = (uint8_t)(data[u] >> 24);
		sha1_update(sc, tmp, 4);
	}
}

static void
sha1_out(const sha1_context *cc, void *dst)
{
	uint8_t tmp[64];
	uint32_t val[5];

	size_t ptr = (size_t)cc->count & 63;
	memcpy(tmp, cc->buf, ptr);
	memcpy(val, cc->val, sizeof val);
	tmp[ptr ++] = 0x80;
	if (ptr > 56) {
		memset(tmp + ptr, 0, 64 - ptr);
		sha1_round_inner(tmp, val);
		memset(tmp, 0, 56);
	} else {
		memset(tmp + ptr, 0, 56 - ptr);
	}
	uint64_t v = cc->count << 3;
	for (int i = 0; i < 8; i ++) {
		tmp[56 + i] = (uint8_t)(v >> (8 * (7 - i)));
	}
	sha1_round_inner(tmp, val);
	uint8_t *buf = dst;
	for (int i = 0; i < 5; i ++) {
		uint32_t w = val[i];
		buf[(i << 2) + 0] = (uint8_t)(w >> 24);
		buf[(i << 2) + 1] = (uint8_t)(w >> 16);
		buf[(i << 2) + 2] = (uint8_t)(w >> 8);
		buf[(i << 2) + 3] = (uint8_t)w;
	}
}

TARGET_AVX2
static void
test_mp31(void)
{
	printf("Test mp31: ");
	fflush(stdout);

	for (unsigned logn = 1; logn <= 10; logn ++) {
		size_t n = (size_t)1 << logn;
		uint32_t *t1 = xmalloc(n * sizeof *t1);
		uint32_t *t2 = xmalloc(n * sizeof *t2);
		uint32_t *w3 = xmalloc(2 * n * sizeof *w3);
		uint32_t *t5 = xmalloc(n * sizeof *t5);

		for (int i = 0; i < 10; i ++) {
			ntrugen_prng_chacha8_context rc;
			uint8_t bb[512];
			uint32_t p = PRIMES[i].p;
			uint32_t p0i = PRIMES[i].p0i;
			uint32_t g = PRIMES[i].g;
			uint32_t ig = PRIMES[i].ig;
			uint32_t R2 = PRIMES[i].R2;

			/* Generate random polynomials p1 and p2. */
			bb[0] = logn;
			bb[1] = i;
			ntrugen_prng_chacha8_init(&rc, bb, 2);
			for (size_t u = 0; u < n; u ++) {
				if ((u & 63) == 0) {
					ntrugen_prng_chacha8_out(&rc,
						bb, sizeof bb);
				}
				size_t v = (u & 63) << 3;
				t1[u] = ((uint32_t)bb[v + 0]
					| ((uint32_t)bb[v + 1] << 8)
					| ((uint32_t)bb[v + 2] << 16)
					| ((uint32_t)bb[v + 3] << 24)) % p;
				t2[u] = ((uint32_t)bb[v + 4]
					| ((uint32_t)bb[v + 5] << 8)
					| ((uint32_t)bb[v + 6] << 16)
					| ((uint32_t)bb[v + 7] << 24)) % p;
			}

			/* Compute product t1*t2 manually (into t5). */
			memset(w3, 0, 2 * n * sizeof *w3);
			for (size_t u = 0; u < n; u ++) {
				for (size_t v = 0; v < n; v ++) {
					uint64_t z = (uint64_t)t1[u]
						* (uint64_t)t2[v]
						+ (uint64_t)w3[u + v];
					w3[u + v] = (uint32_t)(z % p);
				}
			}
			for (size_t u = 0; u < n; u ++) {
				uint32_t x = w3[u];
				uint32_t y = w3[u + n];
				t5[u] = (y > x) ? (x + p) - y : x - y;
			}

			/* Use the NTT to compute the product t1*t2. */
			uint32_t *gm = w3;
			uint32_t *igm = w3 + n;
			mp_mkgmigm(logn, gm, igm, g, ig, p, p0i);

			if (logn >= 3) {
				uint32_t xg = g;
				for (unsigned j = logn; j < 10; j ++) {
					xg = mp_montymul(xg, xg, p, p0i);
				}
				for (unsigned j = 3; j < logn; j ++) {
					xg = mp_montymul(xg, xg, p, p0i);
				}
				for (size_t u = 0; u < n; u += 8) {
					uint32_t xz[8];
					xz[0] = gm[u];
					xz[4] = mp_montymul(xz[0], xg, p, p0i);
					xz[2] = mp_montymul(xz[4], xg, p, p0i);
					xz[6] = mp_montymul(xz[2], xg, p, p0i);
					xz[1] = mp_montymul(xz[6], xg, p, p0i);
					xz[5] = mp_montymul(xz[1], xg, p, p0i);
					xz[3] = mp_montymul(xz[5], xg, p, p0i);
					xz[7] = mp_montymul(xz[3], xg, p, p0i);
					for (size_t v = 1; v < 8; v ++) {
						if (xz[v] != gm[u + v]) {
							fprintf(stderr,
								"ERR gm\n");
							exit(EXIT_FAILURE);
						}
					}
				}
			}

			mp_NTT(logn, t1, gm, p, p0i);
			mp_NTT(logn, t2, gm, p, p0i);

			for (size_t u = 0; u < n; u ++) {
				t1[u] = mp_montymul(R2, mp_montymul(
					t1[u], t2[u], p, p0i), p, p0i);
			}
			mp_iNTT(logn, t1, igm, p, p0i);

			for (size_t u = 0; u < n; u ++) {
				if (t1[u] != t5[u]) {
					fprintf(stderr, "ERR (%zu) %u / %u\n",
						u, t1[u], t5[u]);
					for (size_t v = 0; v < n; v ++) {
						fprintf(stderr, " %08X", t1[v]);
					}
					fprintf(stderr, "\n");
					for (size_t v = 0; v < n; v ++) {
						fprintf(stderr, " %08X", t5[v]);
					}
					fprintf(stderr, "\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		uint32_t p = PRIMES[0].p;
		uint32_t p0i = PRIMES[0].p0i;
		mp_mkgm(logn, t5, PRIMES[0].g, p, p0i);
		for (size_t u = 0; u < n; u ++) {
			t1[u] = 0;
		}
		t1[1] = 1;
		mp_NTT(logn, t1, t5, p, p0i);
		size_t hn = n >> 1;
		for (size_t j = 0; j < hn; j ++) {
			uint32_t r0 = mp_montymul(t5[hn + j], 1, p, p0i);
			uint32_t r1 = p - r0;
			uint32_t v0 = t1[(j << 1) + 0];
			uint32_t v1 = t1[(j << 1) + 1];
			if (v0 != r0 || v1 != r1) {
				fprintf(stderr,
					"ERR NTT(X) (%zu) %u %u / %u %u\n",
					j << 1, v0, v1, r0, r1);
				exit(EXIT_FAILURE);
			}
		}

		xfree(t1);
		xfree(t2);
		xfree(w3);
		xfree(t5);
		printf(".");
		fflush(stdout);
	}

	printf(" ");
	fflush(stdout);
	uint32_t x = 200000;
	uint32_t y = 0;
	uint32_t p = 2147423699;  /* 194813*73*151 */
	for (int i = 0; i < 11023; i ++, x ++, y ++) {
		uint32_t z = mp_div(x, y, p);
		if (z == 0) {
			if ((y % 73) != 0 && (y % 151) != 0) {
				fprintf(stderr, "ERR1: %u / %u -> %u\n",
					x, y, z);
				exit(EXIT_FAILURE);
			}
		} else {
			uint32_t x2 = (uint32_t)(((uint64_t)z * y) % p);
			if (x != x2) {
				fprintf(stderr, "ERR2: %u / %u -> %u\n",
					x, y, z);
				exit(EXIT_FAILURE);
			}
		}
	}
	printf(".");
	fflush(stdout);

	static const int8_t f16_12289_invert[] = {
		-1, -1, 0, -1, 2, -2, 2, 2, 1, -1, -1, 1, 0, 1, 0, -2
	};
	static const int8_t f16_12289_noninvert[] = {
		1, -1, 0, 2, 2, -1, -2, -2, 1, 1, 1, 0, -2, 0, -1, 0
	};
	uint32_t tmp[32];
	if (!poly_is_invertible(4, f16_12289_invert, 2147342993,
		4010727823, 547147950, 12289, 2863078533, 45, tmp))
	{
		fprintf(stderr, "ERR: invertible wrongly dismissed\n");
		exit(EXIT_FAILURE);
	}
	if (poly_is_invertible(4, f16_12289_noninvert, 2147342993,
		4010727823, 547147950, 12289, 2863078533, 45, tmp))
	{
		fprintf(stderr, "ERR: non-invertible wrongly accepted\n");
		exit(EXIT_FAILURE);
	}
	printf(".");
	fflush(stdout);

#if NTRUGEN_AVX2
	ntrugen_prng_chacha8_context rc;
	ntrugen_prng_chacha8_init(&rc, NULL, 0);
	for (int i = 0; i < 30; i ++) {
		p = PRIMES[i].p;
		uint32_t p0i = PRIMES[i].p0i;
		__m256i yp = _mm256_set1_epi32(p);
		__m256i yp0i = _mm256_set1_epi32(p0i);
		uint8_t buf[512];
		uint32_t a[64];
		uint32_t b[64];
		ntrugen_prng_chacha8_out(&rc, buf, sizeof buf);
		for (int j = 0; j < 64; j ++) {
			a[j] = ((uint32_t)buf[(j << 2) + 0]
				| ((uint32_t)buf[(j << 2) + 1] << 8)
				| ((uint32_t)buf[(j << 2) + 1] << 16)
				| ((uint32_t)buf[(j << 2) + 1] << 24)) % p;
			b[j] = ((uint32_t)buf[(j << 2) + 256]
				| ((uint32_t)buf[(j << 2) + 257] << 8)
				| ((uint32_t)buf[(j << 2) + 258] << 16)
				| ((uint32_t)buf[(j << 2) + 259] << 24)) % p;
		}
		for (int j = 0; j < 64; j += 8) {
			union {
				__m256i y;
				uint32_t v[8];
			} z;
			__m256i ya = _mm256_loadu_si256((__m256i *)(a + j));
			__m256i yb = _mm256_loadu_si256((__m256i *)(b + j));
			z.y = mp_montymul_x4(ya, yb, yp, yp0i);
			for (int k = 0; k < 8; k += 2) {
				if (z.v[k] != mp_montymul(
					a[j + k], b[j + k], p, p0i)
					|| z.v[k + 1] != 0)
				{
					fprintf(stderr, "ERR (x4): j=%d k=%d\n",
						j, k);
					fprintf(stderr, "a = 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X\n",
						a[j + 0], a[j + 1],
						a[j + 2], a[j + 3],
						a[j + 4], a[j + 5],
						a[j + 6], a[j + 7]);
					fprintf(stderr, "b = 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X\n",
						b[j + 0], b[j + 1],
						b[j + 2], b[j + 3],
						b[j + 4], b[j + 5],
						b[j + 6], b[j + 7]);
					fprintf(stderr, "c = 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X\n",
						z.v[0], z.v[1],
						z.v[2], z.v[3],
						z.v[4], z.v[5],
						z.v[6], z.v[7]);
					exit(EXIT_FAILURE);
				}
			}
			z.y = mp_montymul_x8(ya, yb, yp, yp0i);
			for (int k = 0; k < 8; k ++) {
				if (z.v[k] != mp_montymul(
					a[j + k], b[j + k], p, p0i))
				{
					fprintf(stderr, "ERR (x8): j=%d k=%d\n",
						j, k);
					fprintf(stderr, "a = 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X\n",
						a[j + 0], a[j + 1],
						a[j + 2], a[j + 3],
						a[j + 4], a[j + 5],
						a[j + 6], a[j + 7]);
					fprintf(stderr, "b = 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X\n",
						b[j + 0], b[j + 1],
						b[j + 2], b[j + 3],
						b[j + 4], b[j + 5],
						b[j + 6], b[j + 7]);
					fprintf(stderr, "c = 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X 0x%08X\n",
						z.v[0], z.v[1],
						z.v[2], z.v[3],
						z.v[4], z.v[5],
						z.v[6], z.v[7]);
					exit(EXIT_FAILURE);
				}
			}
		}
	}
#endif // NTRUGEN_AVX2

	printf(" done.\n");
	fflush(stdout);
}

/*
 * KAT format: each test is: len, x, y, u, v
 * (len is over one word; x, y, u and v are over len words each; terminator
 * is len == 0).
 */
static const uint32_t KAT_BEZOUT[] = {
	1,
	0x7F8EB937,
	0x63FA13AD,
	0x17F61B02,
	0x1E9252C1,
	2,
	0x246DE65D, 0x3D297D77,
	0x4634817B, 0x3504A434,
	0x12CCC029, 0x23C130A3,
	0x1A7B7F6C, 0x293F2E6A,
	3,
	0x5C4094E9, 0x42E26ADC, 0x12F51CB6,
	0x329D6CF1, 0x18F07A5A, 0x0461C01F,
	0x35558751, 0x1A4E04CA, 0x00835687,
	0x55B5A838, 0x18E9084F, 0x023837A2,
	4,
	0x1F8CA55D, 0x5B07B027, 0x0D392DCC, 0x4CAEBF0C,
	0x2DA460E5, 0x7F005900, 0x23C98C1E, 0x231BA1AF,
	0x77A3C826, 0x0872A236, 0x234DF600, 0x1E8C3905,
	0x3275A0C9, 0x62B40B41, 0x0B04EA3E, 0x42B8E59B,
	5,
	0x7901DDAF, 0x41435DC7, 0x31131F1B, 0x7B3F8A25, 0x4268BFEA,
	0x66262653, 0x1FF24232, 0x0523EBE1, 0x0BFEA465, 0x0C109627,
	0x72AA84C4, 0x11004822, 0x2D6A7D14, 0x0D211339, 0x0A6BC3F5,
	0x33CE2CB9, 0x4B9C240E, 0x527B452A, 0x33509F59, 0x395C6422,
	0
};

static void
test_bezout(void)
{
	printf("Test zint_bezout: ");
	fflush(stdout);

	const uint32_t *vp = KAT_BEZOUT;
	for (;;) {
		size_t xlen = *vp ++;
		if (xlen == 0) {
			break;
		}
		const uint32_t *x = vp;
		vp += xlen;
		const uint32_t *y = vp;
		vp += xlen;
		const uint32_t *ru = vp;
		vp += xlen;
		const uint32_t *rv = vp;
		vp += xlen;
		uint32_t u[10], v[10], tmp[40];
		if (zint_bezout(u, v, x, y, xlen, tmp) != 1) {
			fprintf(stderr, "FAIL\n");
			exit(EXIT_FAILURE);
		}
		if (memcmp(u, ru, xlen * sizeof *ru) != 0
			|| memcmp(v, rv, xlen * sizeof *rv) != 0)
		{
			fprintf(stderr, "WRONG:\n");
			fprintf(stderr, "u  =");
			for (size_t j = 0; j < xlen; j ++) {
				fprintf(stderr, " %08X", u[j]);
			}
			fprintf(stderr, "\n");
			fprintf(stderr, "ru =");
			for (size_t j = 0; j < xlen; j ++) {
				fprintf(stderr, " %08X", ru[j]);
			}
			fprintf(stderr, "\n");
			fprintf(stderr, "v  =");
			for (size_t j = 0; j < xlen; j ++) {
				fprintf(stderr, " %08X", v[j]);
			}
			fprintf(stderr, "\n");
			fprintf(stderr, "rv =");
			for (size_t j = 0; j < xlen; j ++) {
				fprintf(stderr, " %08X", rv[j]);
			}
			fprintf(stderr, "\n");
			exit(EXIT_FAILURE);
		}
		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_FFT(void)
{
	printf("Test FFT: ");
	fflush(stdout);

	for (unsigned logn = 1; logn <= 10; logn ++) {
		size_t n = (size_t)1 << logn;
		int8_t *f1 = xmalloc(n);
		int8_t *f2 = xmalloc(n);
		int32_t *f3 = xmalloc(n * sizeof *f3);
		fxr *t1 = xmalloc(n * sizeof *t1);
		fxr *t2 = xmalloc(n * sizeof *t2);

		for (int i = 0; i < 10; i ++) {
			ntrugen_prng_chacha8_context rc;
			uint8_t bb[512];

			/* Generate random polynomials f1 and f2. */
			bb[0] = logn;
			bb[1] = i;
			ntrugen_prng_chacha8_init(&rc, bb, 2);
			ntrugen_prng_chacha8_out(&rc, f1, n);
			ntrugen_prng_chacha8_out(&rc, f2, n);

			/* Test FFT/iFFT round trip on f1. */
			vect_set(logn, t1, f1);
			vect_FFT(logn, t1);
			vect_iFFT(logn, t1);
			for (size_t u = 0; u < n; u ++) {
				if (fxr_round(t1[u]) != f1[u]) {
					fprintf(stderr, "ERR1 (%zu) %d / %d\n",
						u, f1[u], fxr_round(t1[u]));
					exit(EXIT_FAILURE);
				}
			}

			/* Compute product f1*f2 manually (into f3). */
			memset(f3, 0, n * sizeof *f3);
			for (size_t u = 0; u < n; u ++) {
				int32_t x1 = f1[u];
				for (size_t v = 0; v < n; v ++) {
					int32_t x2 = f2[v];
					if (u + v < n) {
						f3[u + v] += x1 * x2;
					} else {
						f3[u + v - n] -= x1 * x2;
					}
				}
			}

			/* Use the FFT to compute the product f1*f2. */
			vect_set(logn, t1, f1);
			vect_set(logn, t2, f2);
			vect_FFT(logn, t1);
			vect_FFT(logn, t2);
			vect_mul_fft(logn, t1, t2);
			vect_iFFT(logn, t1);
			for (size_t u = 0; u < n; u ++) {
				if (f3[u] != fxr_round(t1[u])) {
					fprintf(stderr, "ERR2 (%zu) %d / %d\n",
						u, f3[u], fxr_round(t1[u]));
					exit(EXIT_FAILURE);
				}
			}
		}

		xfree(f1);
		xfree(f2);
		xfree(f3);
		xfree(t1);
		xfree(t2);
		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

/* Defined in ng_bat.h */
int Zn(poly_is_invertible_mod257)(unsigned logn,
	const int8_t *restrict f, uint32_t *restrict tmp);
int Zn(poly_is_invertible_mod769)(unsigned logn,
	const int8_t *restrict f, uint32_t *restrict tmp);

static void
test_inv_mq(void)
{
	printf("Test inv_mq: ");
	fflush(stdout);

	int8_t *f = xmalloc(1024);
	uint16_t *t1 = xmalloc(1025 * sizeof *t1);
	uint16_t *t2 = xmalloc(1025 * sizeof *t2);
	uint32_t *tmp = xmalloc(3072 * sizeof *tmp);
	ntrugen_prng_chacha8_context rng;
	ntrugen_prng_chacha8_init(&rng, "inv_mq", 6);

	for (int i = 0; i < 40; i ++) {
		unsigned logn = (i < 20) ? 9 : 10;
		size_t n = (size_t)1 << logn;
		uint32_t q = (logn == 9) ? 257 : 769;

		if ((i & 1) == 1) {
			uint8_t buf[30];
			memset(f, 0, n);
			ntrugen_prng_chacha8_out(&rng, buf, sizeof buf);
			for (size_t u = 0; u < sizeof buf; u += 2) {
				size_t v = dec16le(buf + u) & (n - 1);
				int nm = 1 - 2 * (buf[u + 1] >> 7);
				if (q == 257) {
					f[v] += 3 * nm;
					if ((v + 4) >= n) {
						f[v + 4 - n] -= nm;
					} else {
						f[v + 4] += nm;
					}
				} else {
					f[v] += 7 * nm;
					if ((v + 8) >= n) {
						f[v + 8 - n] -= nm;
					} else {
						f[v + 8] += nm;
					}
				}
			}
		} else {
			ntrugen_prng_chacha8_out(&rng, f, n);
		}

		/*
		 * A custom GCD over polynomials (mod q).
		 */
		for (size_t u = 0; u < n; u ++) {
			int x = f[u];
			t1[u] = (uint16_t)(x < 0 ? (int)q + x : x);
		}
		memset(t2, 0, n * sizeof *t2);
		t2[0] = 1;
		t2[n] = 1;
		size_t t1len = n;
		while (t1len > 0 && t1[t1len - 1] == 0) {
			t1len --;
		}
		size_t t2len = n + 1;
		while (t1len > 0) {
			if (t1[0] == 0) {
				memmove(t1, t1 + 1, (t1len - 1) * sizeof *t1);
				t1len --;
				continue;
			}
			if (t1len < t2len) {
				for (size_t u = 0; u < t1len; u ++) {
					int x = t1[u];
					t1[u] = t2[u];
					t2[u] = x;
				}
				memmove(t1 + t1len, t2 + t1len,
					(t2len - t1len) * sizeof *t1);
				size_t v = t1len;
				t1len = t2len;
				t2len = v;
			}
			uint32_t m1 = q - t2[0];
			uint32_t m2 = t1[0];
			for (size_t u = 1; u < t2len; u ++) {
				uint32_t x1 = t1[u];
				uint32_t x2 = t2[u];
				t1[u] = (uint16_t)((x1 * m1 + x2 * m2) % q);
			}
			for (size_t u = t2len; u < t1len; u ++) {
				uint32_t x1 = t1[u];
				t1[u] = (uint16_t)((x1 * m1) % q);
			}
			t1[0] = 0;
			while (t1len > 0 && t1[t1len - 1] == 0) {
				t1len --;
			}
		}
		int is_inv_1 = (t2len == 1);

		int is_inv_2 = (q == 257)
			? Zn(poly_is_invertible_mod257)(logn, f, tmp)
			: Zn(poly_is_invertible_mod769)(logn, f, tmp);
		if (is_inv_1 != is_inv_2) {
			fprintf(stderr, "ERR: %d / %d\n", is_inv_1, is_inv_2);
			exit(EXIT_FAILURE);
		}

		printf(".");
		fflush(stdout);
	}

	xfree(f);
	xfree(t1);
	xfree(t2);
	xfree(tmp);

	printf(" done.\n");
	fflush(stdout);
}

static void
test_BAT_keygen(void)
{
	printf("Test BAT_keygen: ");
	fflush(stdout);

	int8_t *f = xmalloc(1024);
	int8_t *g = xmalloc(1024);
	int8_t *F = xmalloc(1024);
	int8_t *G = xmalloc(1024);
	int8_t *G2 = xmalloc(1024);
	int32_t *w = xmalloc(1024 * sizeof *w);
	void *tmpbytes = xmalloc(20 * 1024 + 18);
	uint32_t *tmp = (void *)(((uintptr_t)tmpbytes + 7) & ~(uintptr_t)7);
	sha1_context sc;
	sha1_init(&sc);
	for (unsigned logn = 8; logn <= 10; logn ++) {
		for (int i = 0; i < 10; i ++) {
			ntrugen_prng_chacha8_context rc;
			uint8_t buf[2];
			size_t n = (size_t)1 << logn;
			buf[0] = logn;
			buf[1] = (uint8_t)i;
			ntrugen_prng_chacha8_init(&rc, buf, 2);
			tmp[5 * n] = 0xA7C083FE;
			tmp[5 * n + 1] = 0xA7C083FE;
			tmp[5 * n + 2] = 0xA7C083FE;
			if (i == 0) {
				if (BAT_keygen(logn, f, g, F, NULL, w,
					&ntrugen_prng_chacha8_out, &rc,
					tmp, 20 * n) != 0)
				{
					fprintf(stderr, "ERR: call failed\n");
					exit(EXIT_FAILURE);
				}
				if (tmp[5 * n] != 0xA7C083FE) {
					fprintf(stderr, "ERR: buf overflow\n");
					exit(EXIT_FAILURE);
				}
				if (BAT_recover_G(logn,
					G, f, g, F, tmp, 12 * n) != 0)
				{
					fprintf(stderr, "RECOVER G ERR\n");
					exit(EXIT_FAILURE);
				}

				/* Try again, to test proper generation of w
				   in the buffer */
				buf[0] = logn;
				buf[1] = (uint8_t)i;
				ntrugen_prng_chacha8_init(&rc, buf, 2);
				memset(tmp, 0, 20 * n);
				if (BAT_keygen(logn, f, g, NULL, NULL, NULL,
					&ntrugen_prng_chacha8_out, &rc,
					tmp, 20 * n) != 0)
				{
					fprintf(stderr, "ERR: call 2 failed\n");
					exit(EXIT_FAILURE);
				}
				if (memcmp((int8_t *)tmp + 2 * n, w,
					n * sizeof *w) != 0)
				{
					fprintf(stderr, "ERR: different w\n");
					exit(EXIT_FAILURE);
				}
			} else if (i == 1) {
				if (BAT_keygen(logn, f, g, NULL, NULL, NULL,
					&ntrugen_prng_chacha8_out, &rc,
					(int8_t *)tmp + 1, 20 * n + 7) != 0)
				{
					fprintf(stderr, "ERR: call failed\n");
					exit(EXIT_FAILURE);
				}
				if (tmp[5 * n + 2] != 0xA7C083FE) {
					fprintf(stderr, "ERR: buf overflow\n");
					exit(EXIT_FAILURE);
				}
				memcpy(F, (int8_t *)tmp + 1, n);
				memcpy(G, (int8_t *)tmp + 1 + n, n);
				if (BAT_recover_G(logn,
					NULL, f, g, F,
					(int8_t *)tmp + 1, 12 * n + 7) != 0)
				{
					fprintf(stderr, "RECOVER G ERR\n");
					exit(EXIT_FAILURE);
				}
				if (memcmp(G, (int8_t *)tmp + 1, n) != 0) {
					fprintf(stderr, "WRONG RECOVERED G\n");
					exit(EXIT_FAILURE);
				}
				tmp[5 * n] = 0xA7C083FE;
				tmp[5 * n + 1] = 0xA7C083FE;
			} else {
				if (BAT_keygen(logn, f, g, F, G2, w,
					&ntrugen_prng_chacha8_out, &rc,
					tmp, 20 * n) != 0)
				{
					fprintf(stderr, "ERR: call failed\n");
					exit(EXIT_FAILURE);
				}
				if (tmp[5 * n] != 0xA7C083FE) {
					fprintf(stderr, "ERR: buf overflow\n");
					exit(EXIT_FAILURE);
				}
				if (BAT_recover_G(logn,
					G, f, g, F, tmp, 12 * n) != 0)
				{
					fprintf(stderr, "RECOVER G ERR\n");
					exit(EXIT_FAILURE);
				}
				if (memcmp(G, G2, n) != 0) {
					fprintf(stderr, "WRONG RECOVERED G\n");
					exit(EXIT_FAILURE);
				}
			}
			for (size_t u = 0; u < n; u ++) {
				int32_t s = 0;
				for (size_t v = 0; v <= u; v ++) {
					s += (int32_t)f[v] * G[u - v];
					s -= (int32_t)g[v] * F[u - v];
				}
				for (size_t v = u + 1; v < n; v ++) {
					s -= (int32_t)f[v] * G[n + u - v];
					s += (int32_t)g[v] * F[n + u - v];
				}
				int32_t rv = 0;
				if (u == 0) {
					switch (logn) {
					case 8: rv = -128; break;
					case 9: rv = -257; break;
					case 10: rv = -769; break;
					}
				}
				if (rv != s) {
					fprintf(stderr, "KEYGEN ERR\n");
					exit(EXIT_FAILURE);
				}
			}

			sha1_update(&sc, f, n);
			sha1_update(&sc, g, n);
			sha1_update(&sc, F, n);
			sha1_update(&sc, G, n);
		}

		printf(".");
		fflush(stdout);
	}
	xfree(f);
	xfree(g);
	xfree(F);
	xfree(G);
	xfree(G2);
	xfree(w);
	xfree(tmpbytes);

	uint8_t hv[20], hvref[20];
	sha1_out(&sc, hv);
	HEXTOBIN(hvref, "3a71dc372bb73a6f4b86eeea5dc2b553f2afeec4");
	check_eq(hv, hvref, sizeof hv, "aggregate hash");

	printf(" done.\n");
	fflush(stdout);
}

static void
test_Falcon_keygen(void)
{
	printf("Test Falcon_keygen: ");
	fflush(stdout);

	int8_t *f = xmalloc(1024);
	int8_t *g = xmalloc(1024);
	int8_t *F = xmalloc(1024);
	int8_t *G = xmalloc(1024);
	void *tmpbytes = xmalloc(20 * 1024 + 18);
	uint32_t *tmp = (void *)(((uintptr_t)tmpbytes + 7) & ~(uintptr_t)7);
	sha1_context sc;
	sha1_init(&sc);
	for (unsigned logn = 2; logn <= 10; logn ++) {
		for (int i = 0; i < 10; i ++) {
			ntrugen_prng_chacha8_context rc;
			uint8_t buf[2];
			size_t n = (size_t)1 << logn;
			buf[0] = logn;
			buf[1] = (uint8_t)i;
			ntrugen_prng_chacha8_init(&rc, buf, 2);
			tmp[5 * n] = 0xA7C083FE;
			tmp[5 * n + 1] = 0xA7C083FE;
			tmp[5 * n + 2] = 0xA7C083FE;
			if (i == 0) {
				if (Falcon_keygen(logn, f, g, F, NULL,
					&ntrugen_prng_chacha8_out, &rc,
					tmp, 20 * n) != 0)
				{
					fprintf(stderr, "ERR: call failed\n");
					exit(EXIT_FAILURE);
				}
				if (tmp[5 * n] != 0xA7C083FE) {
					fprintf(stderr, "ERR: buf overflow\n");
					exit(EXIT_FAILURE);
				}
				memcpy(G, (int8_t *)tmp + n, n);
			} else if (i == 1) {
				if (Falcon_keygen(logn, f, g, NULL, NULL,
					&ntrugen_prng_chacha8_out, &rc,
					(int8_t *)tmp + 1, 20 * n + 7) != 0)
				{
					fprintf(stderr, "ERR: call failed\n");
					exit(EXIT_FAILURE);
				}
				if (tmp[5 * n + 2] != 0xA7C083FE) {
					fprintf(stderr, "ERR: buf overflow\n");
					exit(EXIT_FAILURE);
				}
				memcpy(F, (int8_t *)tmp + 1, n);
				memcpy(G, (int8_t *)tmp + 1 + n, n);
				tmp[5 * n] = 0xA7C083FE;
				tmp[5 * n + 1] = 0xA7C083FE;
			} else {
				if (Falcon_keygen(logn, f, g, F, G,
					&ntrugen_prng_chacha8_out, &rc,
					tmp, 20 * n) != 0)
				{
					fprintf(stderr, "ERR: call failed\n");
					exit(EXIT_FAILURE);
				}
				if (tmp[5 * n] != 0xA7C083FE) {
					fprintf(stderr, "ERR: buf overflow\n");
					exit(EXIT_FAILURE);
				}
			}
			for (size_t u = 0; u < n; u ++) {
				int32_t s = 0;
				for (size_t v = 0; v <= u; v ++) {
					s += (int32_t)f[v] * G[u - v];
					s -= (int32_t)g[v] * F[u - v];
				}
				for (size_t v = u + 1; v < n; v ++) {
					s -= (int32_t)f[v] * G[n + u - v];
					s += (int32_t)g[v] * F[n + u - v];
				}
				int32_t rv = (u == 0) ? 12289 : 0;
				if (rv != s) {
					fprintf(stderr, "KEYGEN ERR\n");
					exit(EXIT_FAILURE);
				}
			}

			sha1_update(&sc, f, n);
			sha1_update(&sc, g, n);
			sha1_update(&sc, F, n);
			sha1_update(&sc, G, n);
		}

		printf(".");
		fflush(stdout);
	}
	xfree(f);
	xfree(g);
	xfree(F);
	xfree(G);
	xfree(tmpbytes);

	uint8_t hv[20], hvref[20];
	sha1_out(&sc, hv);
	HEXTOBIN(hvref, "9e07610b81a57c2eb8a7384ad29dd7b10d00c853");
	check_eq(hv, hvref, sizeof hv, "aggregate hash");

	printf(" done.\n");
	fflush(stdout);
}

static void
verify_pub_Hawk(unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const int16_t *restrict q00, const int16_t *restrict q01,
	const int32_t *restrict q11)
{
	/*
	 * Check range.
	 */
	size_t n = (size_t)1 << logn;
	int lim00, lim01;
	int32_t lim11;
	if (logn == 8) {
		lim00 = 1 << 9;
		lim01 = 1 << 11;
		lim11 = (int32_t)1 << 13;
	} else if (logn == 9) {
		lim00 = 1 << 9;
		lim01 = 1 << 12;
		lim11 = (int32_t)1 << 15;
	} else {
		lim00 = 1 << 10;
		lim01 = 1 << 14;
		lim11 = (int32_t)1 << 17;
	}
	for (size_t u = 1; u < n; u ++) {
		if (q00[u] <= -lim00 || q00[u] >= +lim00) {
			fprintf(stderr,
				"ERR RANGE: q00[%zu] = %d\n", u, q00[u]);
			exit(EXIT_FAILURE);
		}
		if (q11[u] <= -lim11 || q11[u] >= +lim11) {
			fprintf(stderr,
				"ERR RANGE: q11[%zu] = %d\n", u, q11[u]);
			exit(EXIT_FAILURE);
		}
	}
	for (size_t u = 0; u < n; u ++) {
		if (q01[u] <= -lim01 || q01[u] >= +lim01) {
			fprintf(stderr,
				"ERR RANGE: q01[%zu] = %d\n", u, q01[u]);
			exit(EXIT_FAILURE);
		}
	}

	/*
	 * q00 = f*adj(f) + g*adj(g)
	 */
	for (size_t u = 0; u < n; u ++) {
		int32_t s = (int32_t)f[0] * (int32_t)f[u];
		s += (int32_t)g[0] * (int32_t)g[u];
		for (size_t v = 1; v <= u; v ++) {
			s -= (int32_t)f[n - v] * (int32_t)f[u - v];
			s -= (int32_t)g[n - v] * (int32_t)g[u - v];
		}
		for (size_t v = u + 1; v < n; v ++) {
			s += (int32_t)f[n - v] * (int32_t)f[n + u - v];
			s += (int32_t)g[n - v] * (int32_t)g[n + u - v];
		}
		if (s != q00[u]) {
			fprintf(stderr, "ERR q00[%zu] -> %d / %d\n",
				u, s, q00[u]);
			exit(EXIT_FAILURE);
		}
	}

	/*
	 * q01 = F*adj(f) + G*adj(g)
	 */
	for (size_t u = 0; u < n; u ++) {
		int32_t s = (int32_t)f[0] * (int32_t)F[u];
		s += (int32_t)g[0] * (int32_t)G[u];
		for (size_t v = 1; v <= u; v ++) {
			s -= (int32_t)f[n - v] * (int32_t)F[u - v];
			s -= (int32_t)g[n - v] * (int32_t)G[u - v];
		}
		for (size_t v = u + 1; v < n; v ++) {
			s += (int32_t)f[n - v] * (int32_t)F[n + u - v];
			s += (int32_t)g[n - v] * (int32_t)G[n + u - v];
		}
		if (s != q01[u]) {
			fprintf(stderr, "ERR q01[%zu] -> %d / %d\n",
				u, s, q01[u]);
			exit(EXIT_FAILURE);
		}
	}

	/*
	 * q11 = F*adj(F) + G*adj(G)
	 */
	for (size_t u = 0; u < n; u ++) {
		int32_t s = (int32_t)F[0] * (int32_t)F[u];
		s += (int32_t)G[0] * (int32_t)G[u];
		for (size_t v = 1; v <= u; v ++) {
			s -= (int32_t)F[n - v] * (int32_t)F[u - v];
			s -= (int32_t)G[n - v] * (int32_t)G[u - v];
		}
		for (size_t v = u + 1; v < n; v ++) {
			s += (int32_t)F[n - v] * (int32_t)F[n + u - v];
			s += (int32_t)G[n - v] * (int32_t)G[n + u - v];
		}
		if (s != q11[u]) {
			fprintf(stderr, "ERR q11[%zu] -> %d / %d\n",
				u, s, q11[u]);
			exit(EXIT_FAILURE);
		}
	}
}

static void
test_Hawk_keygen(void)
{
	printf("Test Hawk_keygen: ");
	fflush(stdout);

	int8_t *f = xmalloc(1024);
	int8_t *g = xmalloc(1024);
	int8_t *F = xmalloc(1024);
	int8_t *G = xmalloc(1024);
	int8_t *G2 = xmalloc(1024);
	int16_t *q00 = xmalloc(2048);
	int16_t *q01 = xmalloc(2048);
	int32_t *q11 = xmalloc(4096);
	size_t max_epk_len = 3000;
	uint8_t *epk = xmalloc(max_epk_len);
	void *tmpbytes = xmalloc(20 * 1024 + 18);
	uint32_t *tmp = (void *)(((uintptr_t)tmpbytes + 7) & ~(uintptr_t)7);
	sha1_context sc;
	sha1_init(&sc);

	for (unsigned logn = 8; logn <= 10; logn ++) {
		for (int i = 0; i < 10; i ++) {
			shake_context rc;
			uint8_t buf[2];
			size_t n = (size_t)1 << logn;
			buf[0] = logn;
			buf[1] = (uint8_t)i;
			shake_init(&rc, 256);
			shake_inject(&rc, buf, 2);
			shake_flip(&rc);
			tmp[5 * n] = 0xA7C083FE;
			tmp[5 * n + 1] = 0xA7C083FE;
			tmp[5 * n + 2] = 0xA7C083FE;
			uint8_t seed[40];
			memset(seed, (uint8_t)i, sizeof seed);
			size_t seed_len = 8 + ((size_t)1 << (logn - 5));
			if (i == 0) {
				if (Hawk_keygen(logn, f, g, F, NULL,
					q00, q01, q11, seed,
					(ntrugen_rng)&shake_extract, &rc,
					tmp, 20 * n) != 0)
				{
					fprintf(stderr, "ERR: call failed\n");
					exit(EXIT_FAILURE);
				}
				if (tmp[5 * n] != 0xA7C083FE) {
					fprintf(stderr, "ERR: buf overflow\n");
					exit(EXIT_FAILURE);
				}
				if (Hawk_recover_G(logn,
					G, f, g, F, tmp, 12 * n) != 0)
				{
					fprintf(stderr, "RECOVER G ERR\n");
					exit(EXIT_FAILURE);
				}
			} else if (i == 1) {
				if (Hawk_keygen(logn, f, g, NULL, NULL,
					NULL, NULL, NULL, NULL,
					(ntrugen_rng)&shake_extract, &rc,
					(int8_t *)tmp + 1, 20 * n + 7) != 0)
				{
					fprintf(stderr, "ERR: call failed\n");
					exit(EXIT_FAILURE);
				}
				if (tmp[5 * n + 2] != 0xA7C083FE) {
					fprintf(stderr, "ERR: buf overflow\n");
					exit(EXIT_FAILURE);
				}
				memcpy(F, (int8_t *)tmp + 1, n);
				memcpy(G, (int8_t *)tmp + 1 + n, n);
				memcpy(q00, (int16_t *)
					((int8_t *)tmp + 1 + 2 * n), 2 * n);
				memcpy(q01,
					((int8_t *)tmp + 1 + 4 * n), 2 * n);
				memcpy(q11,
					((int8_t *)tmp + 1 + 6 * n), 4 * n);
				memcpy(seed,
					((int8_t *)tmp + 1 + 10 * n), seed_len);
				if (Hawk_recover_G(logn,
					NULL, f, g, F,
					(int8_t *)tmp + 1, 12 * n + 7) != 0)
				{
					fprintf(stderr, "RECOVER G ERR\n");
					exit(EXIT_FAILURE);
				}
				if (memcmp(G, (int8_t *)tmp + 1, n) != 0) {
					fprintf(stderr, "WRONG RECOVERED G\n");
					exit(EXIT_FAILURE);
				}
				tmp[5 * n] = 0xA7C083FE;
				tmp[5 * n + 1] = 0xA7C083FE;
			} else {
				if (Hawk_keygen(logn, f, g, F, G2,
					q00, q01, q11, seed,
					(ntrugen_rng)&shake_extract, &rc,
					tmp, 20 * n) != 0)
				{
					fprintf(stderr, "ERR: call failed\n");
					exit(EXIT_FAILURE);
				}
				if (tmp[5 * n] != 0xA7C083FE) {
					fprintf(stderr, "ERR: buf overflow\n");
					exit(EXIT_FAILURE);
				}
				if (Hawk_recover_G(logn,
					G, f, g, F, tmp, 12 * n) != 0)
				{
					fprintf(stderr, "RECOVER G ERR\n");
					exit(EXIT_FAILURE);
				}
				if (memcmp(G, G2, n) != 0) {
					fprintf(stderr, "WRONG RECOVERED G\n");
					exit(EXIT_FAILURE);
				}
			}
			int8_t *tf = (int8_t *)tmp;
			int8_t *tg = tf + n;
			Hawk_regen_fg(logn, tf, tg, seed);
			check_eq(f, tf, n, "regen f");
			check_eq(g, tg, n, "regen g");
			for (size_t u = 0; u < n; u ++) {
				int32_t s = 0;
				for (size_t v = 0; v <= u; v ++) {
					s += (int32_t)f[v] * G[u - v];
					s -= (int32_t)g[v] * F[u - v];
				}
				for (size_t v = u + 1; v < n; v ++) {
					s -= (int32_t)f[v] * G[n + u - v];
					s += (int32_t)g[v] * F[n + u - v];
				}
				int32_t rv = (u == 0) ? 1 : 0;
				if (rv != s) {
					fprintf(stderr, "KEYGEN ERR\n");
					exit(EXIT_FAILURE);
				}
			}
			verify_pub_Hawk(logn, f, g, F, G, q00, q01, q11);

			sha1_update(&sc, f, n);
			sha1_update(&sc, g, n);
			sha1_update(&sc, F, n);
			sha1_update(&sc, G, n);
			sha1_update_u16le(&sc, (uint16_t *)q00, n);
			sha1_update_u16le(&sc, (uint16_t *)q01, n);
			sha1_update_u32le(&sc, (uint32_t *)q11, n);

			if (Hawk_recover_qq(logn, NULL, NULL, NULL,
				f, g, F, G, tmp, 18 * n) != 0)
			{
				fprintf(stderr, "RECOVER QQ ERR\n");
				exit(EXIT_FAILURE);
			}
			if (memcmp(q00, tmp, n * 2) != 0) {
				fprintf(stderr, "WRONG RECOVERED q00\n");
				exit(EXIT_FAILURE);
			}
			if (memcmp(q01, (uint8_t *)tmp + 2 * n, n * 2) != 0) {
				fprintf(stderr, "WRONG RECOVERED q01\n");
				exit(EXIT_FAILURE);
			}
			if (memcmp(q11, (uint8_t *)tmp + 4 * n, n * 4) != 0) {
				fprintf(stderr, "WRONG RECOVERED q11\n");
				exit(EXIT_FAILURE);
			}
		}

		printf(".");
		fflush(stdout);
	}
	xfree(f);
	xfree(g);
	xfree(F);
	xfree(G);
	xfree(G2);
	xfree(q00);
	xfree(q01);
	xfree(q11);
	xfree(epk);
	xfree(tmpbytes);

	uint8_t hv[20], hvref[20];
	sha1_out(&sc, hv);
	HEXTOBIN(hvref, "8efd06eb541a2192fa4c407ecb87cf9690204cec");
	check_eq(hv, hvref, sizeof hv, "aggregate hash");

	printf(" done.\n");
	fflush(stdout);
}

#if NTRUGEN_STATS
static int
eqstr_nocase(const char *s1, const char *s2)
{
	for (;;) {
		int c1 = *s1 ++;
		int c2 = *s2 ++;
		if (c1 == 0 || c2 == 0) {
			return c1 == 0 && c2 == 0;
		}
		if (c1 >= 'A' && c1 <= 'Z') {
			c1 += 'a' - 'A';
		}
		if (c2 >= 'A' && c2 <= 'Z') {
			c2 += 'a' - 'A';
		}
		if (c1 != c2) {
			return 0;
		}
	}
}
#endif

int
main(int argc, char *argv[])
{
	(void)argc;
	(void)argv;
#if NTRUGEN_STATS
	stats_init();
#endif
	test_mp31();
	test_bezout();
	test_FFT();
	test_inv_mq();
	test_BAT_keygen();
	test_Falcon_keygen();
	test_Hawk_keygen();

#if NTRUGEN_STATS
	if (argc > 1) {
		unsigned logn;
		int algo;
		const char *algo_name;
		if (argc == 2) {
			logn = strtoul(argv[1], 0, 0);
			algo = 2;
			algo_name = "Falcon";
		} else {
			if (eqstr_nocase(argv[1], "BAT")) {
				algo = 1;
				algo_name = "BAT";
			} else if (eqstr_nocase(argv[1], "Falcon")) {
				algo = 2;
				algo_name = "Falcon";
			} else if (eqstr_nocase(argv[1], "Hawk")) {
				algo = 3;
				algo_name = "Hawk";
			} else {
				fprintf(stderr, "unknown algorithm name\n");
				exit(EXIT_FAILURE);
			}
			logn = strtoul(argv[2], 0, 0);
		}
		if (logn < 2 || logn > 10) {
			fprintf(stderr, "bad log(degree)\n");
			exit(EXIT_FAILURE);
		}
		printf("%s, n=%zu:\n", algo_name, (size_t)1 << logn);
		stats_init();
		size_t n = (size_t)1 << logn;
		int8_t *f = xmalloc(n);
		int8_t *g = xmalloc(n);
		int8_t *F = xmalloc(n);
		int8_t *G = xmalloc(n);
		int32_t *w = (algo == 1) ? xmalloc(n * sizeof *w) : NULL;
		size_t tmp_len = ((size_t)24 << logn) + 7;
		void *tmp = xmalloc(tmp_len);
		ntrugen_prng_chacha8_context rc;
		ntrugen_prng_chacha8_init(&rc, NULL, 0);
		for (uint32_t j = 1; j <= 10000; j ++) {
			int status;
			switch (algo) {
			case 1:
				status = BAT_keygen(logn, f, g, F, G, w,
					&ntrugen_prng_chacha8_out, &rc,
					tmp, tmp_len);
				break;
			case 2:
				status = Falcon_keygen(logn, f, g, F, G,
					&ntrugen_prng_chacha8_out, &rc,
					tmp, tmp_len);
				break;
			case 3:
				status = Hawk_keygen(logn, f, g, F, G,
					NULL, NULL, NULL, NULL,
					&ntrugen_prng_chacha8_out, &rc,
					tmp, tmp_len);
				break;
			default:
				exit(EXIT_FAILURE);
			}
			if (status != 0) {
				fprintf(stderr, "ntrugen failed: %d\n", status);
				exit(EXIT_FAILURE);
			}
			if (j % 100 == 0) {
				uint32_t sa = stats_solve_attempt;
				uint32_t seg = stats_solve_err_gcd;
				uint32_t ser = stats_solve_err_reduce;
				uint32_t sel = stats_solve_err_limit;
				uint32_t ss = stats_solve_success;
				printf("%5u: %5u  ok=%.3f  solve=%.3f"
					" (g=%.3f r=%.3f l=%.3f)",
					j,
					sa,
					(j * 100.0) / sa,
					(ss * 100.0) / sa,
					(seg * 100.0) / sa,
					(ser * 100.0) / sa,
					(sel * 100.0) / sa);
				if (algo == 3) {
					uint32_t ca = stats_hawk_ctt_attempt;
					uint32_t cr = stats_hawk_ctt_reject;
					printf("  cttrej=%.3f",
						(cr * 100.0) / ca);
				}
				printf("\n");
				if (algo == 1) {
					uint32_t wa = stats_compute_w_attempt;
					uint32_t we1 = stats_compute_w_err_lim1;
					uint32_t we2 = stats_compute_w_err_lim2;
					uint32_t we3 = stats_compute_w_err_lim3;
					uint32_t wen = stats_compute_w_err_norm;
					uint32_t ws = stats_compute_w_success;
					printf("                             "
						" w: %.3f (l1=%.3f l2=%.3f"
						" l3=%.3f dn=%.3f)\n",
						(ws * 100.0) / wa,
						(we1 * 100.0) / wa,
						(we2 * 100.0) / wa,
						(we3 * 100.0) / wa,
						(wen * 100.0) / wa);
				}
			}
		}
		xfree(f);
		xfree(g);
		xfree(F);
		xfree(G);
		xfree(w);
		xfree(tmp);
	}
#endif

	return 0;
}
