#include "ng_inner.h"

/*
 * Table format: first entry is kmax; it is followed by 2*kmax values
 * (formal table entries for -kmax to +kmax-1).
 */

/*
 * Type for a PRNG. When invoked, it writes len bytes into the buffer
 * pointed to by dst; the context parameter points to the PRNG state
 * structure. This is a wrapper around the externally provided RNG,
 * so that the external RNG is invoked only for producing chunks of
 * 512 bytes.
 */
typedef struct {
	uint8_t buf[512];
	size_t ptr;
	ntrugen_rng rng;
	void *rng_context;
} prng_buffer;

static inline void
prng_buffer_init(prng_buffer *pb, ntrugen_rng rng, void *rng_context)
{
	pb->ptr = sizeof pb->buf;
	pb->rng = rng;
	pb->rng_context = rng_context;
}

static inline uint16_t
prng_buffer_next_u16(prng_buffer *pb)
{
	if (pb->ptr > (sizeof pb->buf) - 2) {
		pb->rng(pb->rng_context, pb->buf, sizeof pb->buf);
		pb->ptr = 0;
	}
	unsigned x = pb->buf[pb->ptr ++];
	x |= (unsigned)pb->buf[pb->ptr ++] << 8;
	return x;
}

#if NTRUGEN_AVX2
TARGET_AVX2
static inline __m256i
prng_buffer_next_u16_x8(prng_buffer *pb)
{
	if (pb->ptr > (sizeof pb->buf) - 16) {
		pb->rng(pb->rng_context, pb->buf, sizeof pb->buf);
		pb->ptr = 0;
	}
	__m128i x = _mm_loadu_si128((const __m128i *)(pb->buf + pb->ptr));
	pb->ptr += 16;
	return _mm256_cvtepu16_epi32(x);
}
#endif // NTRUGEN_AVX2

/* see ng_inner.h */
TARGET_AVX2
void
gauss_sample_poly(unsigned logn, int8_t *f,
	const uint16_t *tab, ntrugen_rng rng, void *rng_context)
{
	size_t n = (size_t)1 << logn;
	size_t kmax = tab[0];
	prng_buffer pb;
	prng_buffer_init(&pb, rng, rng_context);
#if NTRUGEN_AVX2
	if (logn >= 3 && kmax <= 17) {
		/* With 34 entries, we cover all BAT and Hawk variants,
		   as well as Falcon-512 and Falcon-1024. */
		__m256i ytab[34];
		for (size_t u = 0; u < (kmax << 1); u ++) {
			ytab[u] = _mm256_set1_epi32(tab[u + 1]);
		}
		__m256i yb = _mm256_set1_epi32(-(int32_t)kmax);
		__m256i ys1 = _mm256_setr_epi8(
			0, 4, 8, 12, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			0, 4, 8, 12, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1);
		__m256i ys2 = _mm256_setr_epi32(
			0, 4, 1, 1, 5, 5, 5, 5);
		for (;;) {
			__m256i ypar = _mm256_setzero_si256();
			for (size_t j = 0; j < n; j += 8) {
				__m256i yv = yb;
				__m256i yr = prng_buffer_next_u16_x8(&pb);
				for (size_t k = 0; k < (kmax << 1); k ++) {
					yv = _mm256_sub_epi32(yv,
						_mm256_cmpgt_epi32(
							yr, ytab[k]));
				}
				ypar = _mm256_xor_si256(ypar, yv);
				yv = _mm256_shuffle_epi8(yv, ys1);
				yv = _mm256_permutevar8x32_epi32(yv, ys2);
				__m128i xv = _mm256_castsi256_si128(yv);
#if defined __x86_64__ || defined _M_X64
				*(uint64_t *)(f + j) = _mm_cvtsi128_si64(xv);
#else
				*(uint32_t *)(f + j) = _mm_cvtsi128_si32(xv);
				*(uint32_t *)(f + j + 4) = _mm_cvtsi128_si32(
					_mm_srli_epi64(xv, 32));
#endif
			}
			uint32_t pm = _mm256_movemask_epi8(
				_mm256_slli_epi32(ypar, 31));
			pm ^= (pm >> 16);
			pm ^= (pm >> 8);
			pm ^= (pm >> 4);
			if ((pm & 0x08) == 0x08) {
				return;
			}
		}
	}
#endif // NTRUGEN_AVX2
	for (;;) {
		uint32_t parity = 0;
		for (size_t j = 0; j < n; j ++) {
			uint32_t v = -(uint32_t)kmax;
			uint32_t x = prng_buffer_next_u16(&pb);
			for (size_t k = 1; k <= (kmax << 1); k ++) {
				v += ((uint32_t)tab[k] - x) >> 31;
			}
			f[j] = (int8_t)*(int32_t *)&v;
			parity ^= v;
		}
		if ((parity & 1) != 0) {
			return;
		}
	}
}

/* see ng_inner.h */
void
gauss_sample_poly_reduced(unsigned logn, int8_t *f,
	const uint16_t *tab, ntrugen_rng rng, void *rng_context)
{
	size_t n = (size_t)1 << logn;
	int g = 1 << (8 - logn);
	size_t kmax = tab[0];
	prng_buffer pb;
	prng_buffer_init(&pb, rng, rng_context);
	for (;;) {
		uint32_t parity = 0;
		for (size_t j = 0; j < n;) {
			uint32_t v = -(uint32_t)kmax << (8 - logn);
			for (int i = 0; i < g; i ++) {
				uint32_t x = prng_buffer_next_u16(&pb);
				for (size_t k = 1; k <= (kmax << 1); k ++) {
					uint32_t z = tab[k];
					v += (z - x) >> 31;
				}
			}
			int32_t y = *(int32_t *)&v;
			if (y < -127 || y > +127) {
				continue;
			}
			f[j ++] = (int8_t)y;
			parity ^= v;
		}
		if ((parity & 1) != 0) {
			return;
		}
	}
}
