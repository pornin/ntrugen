#include "ng_inner.h"

/*
 * Convert source f and g into RNS+NTT, at the start of the provided tmp[]
 * (one word per coefficient).
 *
 * RAM USAGE: 3*(2^logn_top)
 */
static void
make_fg_zero(unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g,
	uint32_t *restrict tmp)
{
	size_t n = (size_t)1 << logn;
	uint32_t *ft = tmp;
	uint32_t *gt = ft + n;
	uint32_t *gm = gt + n;
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	poly_mp_set_small(logn, ft, f, p);
	poly_mp_set_small(logn, gt, g, p);
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	mp_NTT(logn, ft, gm, p, p0i);
	mp_NTT(logn, gt, gm, p, p0i);
}

/*
 * One step of computing (f,g) at a given depth.
 * Input: (f,g) of degree 2^(logn_top-depth)
 * Output: (f',g') of degree 2^(logn_top-(depth+1))
 * Input and output values are at the start of tmp[], in RNS+NTT notation.
 *
 * RAM USAGE: 3*(2^logn_top) (at most)
 * (assumptions: max_bl_small[0] = max_bl_small[1] = 1, max_bl_small[2] = 2)
 */
TARGET_AVX2
static void
make_fg_step(const ntru_profile *prof,
	unsigned logn_top, unsigned depth, uint32_t *tmp)
{
	unsigned logn = logn_top - depth;
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t slen = prof->max_bl_small[depth];
	size_t tlen = prof->max_bl_small[depth + 1];

	/*
	 * Layout:
	 *   fd    output f' (hn*tlen)
	 *   gd    output g' (hn*tlen)
	 *   fs    source (n*slen)
	 *   gs    source (n*slen)
	 *   t1    NTT support (n)
	 *   t2    extra (max(n, slen - n))
	 */
	uint32_t *fd = tmp;
	uint32_t *gd = fd + hn * tlen;
	uint32_t *fs = gd + hn * tlen;
	uint32_t *gs = fs + n * slen;
	uint32_t *t1 = gs + n * slen;
	uint32_t *t2 = t1 + n;
	memmove(fs, tmp, 2 * n * slen * sizeof *tmp);

	/*
	 * First slen words: we use the input values directly, and apply
	 * inverse NTT as we go, so that we get the sources in RNS (non-NTT).
	 */
	uint32_t *xf = fs;
	uint32_t *xg = gs;
	uint32_t *yf = fd;
	uint32_t *yg = gd;
	for (size_t u = 0; u < slen; u ++) {
		uint32_t p = PRIMES[u].p;
		uint32_t p0i = PRIMES[u].p0i;
		uint32_t R2 = PRIMES[u].R2;
		for (size_t v = 0; v < hn; v ++) {
			yf[v] = mp_montymul(
				mp_montymul(xf[2 * v], xf[2 * v + 1], p, p0i),
				R2, p, p0i);
			yg[v] = mp_montymul(
				mp_montymul(xg[2 * v], xg[2 * v + 1], p, p0i),
				R2, p, p0i);
		}
		mp_mkigm(logn, t1, PRIMES[u].ig, p, p0i);
		mp_iNTT(logn, xf, t1, p, p0i);
		mp_iNTT(logn, xg, t1, p, p0i);
		xf += n;
		xg += n;
		yf += hn;
		yg += hn;
	}

	/*
	 * Now that fs and gs are in RNS, rebuild their plain integer
	 * coefficients.
	 */
	zint_rebuild_CRT(fs, slen, n, 2, 1, t1);

	/*
	 * Remaining output words.
	 */
	for (size_t u = slen; u < tlen; u ++) {
		uint32_t p = PRIMES[u].p;
		uint32_t p0i = PRIMES[u].p0i;
		uint32_t R2 = PRIMES[u].R2;
		uint32_t Rx = mp_Rx31(slen, p, p0i, R2);
		mp_mkgm(logn, t1, PRIMES[u].g, p, p0i);
#if NTRUGEN_AVX2
		if (logn >= 3) {
			__m256i yp = _mm256_set1_epi32(p);
			__m256i yp0i = _mm256_set1_epi32(p0i);
			__m256i yR2 = _mm256_set1_epi32(R2);
			__m256i yRx = _mm256_set1_epi32(Rx);
			for (size_t v = 0; v < n; v += 8) {
				__m256i yt = zint_mod_small_signed_x8(
					fs + v, slen, n, yp, yp0i, yR2, yRx);
				_mm256_storeu_si256((__m256i *)(t2 + v), yt);
			}
			mp_NTT(logn, t2, t1, p, p0i);
			for (size_t v = 0; v < hn; v += 4) {
				__m256i yt = _mm256_loadu_si256(
					(__m256i *)(t2 + (2 * v)));
				yt = mp_montymul_x4(yt,
					_mm256_srli_epi64(yt, 32), yp, yp0i);
				yt = mp_montymul_x4(yt, yR2, yp, yp0i);
				yt = _mm256_shuffle_epi32(yt, 0xD8);
				yt = _mm256_permute4x64_epi64(yt, 0xD8);
				_mm_storeu_si128((__m128i *)(yf + v),
					_mm256_castsi256_si128(yt));
			}
			yf += hn;
			for (size_t v = 0; v < n; v += 8) {
				__m256i yt = zint_mod_small_signed_x8(
					gs + v, slen, n, yp, yp0i, yR2, yRx);
				_mm256_storeu_si256((__m256i *)(t2 + v), yt);
			}
			mp_NTT(logn, t2, t1, p, p0i);
			for (size_t v = 0; v < hn; v += 4) {
				__m256i yt = _mm256_loadu_si256(
					(__m256i *)(t2 + (2 * v)));
				yt = mp_montymul_x4(yt,
					_mm256_srli_epi64(yt, 32), yp, yp0i);
				yt = mp_montymul_x4(yt, yR2, yp, yp0i);
				yt = _mm256_shuffle_epi32(yt, 0xD8);
				yt = _mm256_permute4x64_epi64(yt, 0xD8);
				_mm_storeu_si128((__m128i *)(yg + v),
					_mm256_castsi256_si128(yt));
			}
			yg += hn;
			continue;
		}
#endif // NTRUGEN_AVX2
		for (size_t v = 0; v < n; v ++) {
			t2[v] = zint_mod_small_signed(
				fs + v, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t2, t1, p, p0i);
		for (size_t v = 0; v < hn; v ++) {
			yf[v] = mp_montymul(
				mp_montymul(t2[2 * v], t2[2 * v + 1], p, p0i),
				R2, p, p0i);
		}
		yf += hn;
		for (size_t v = 0; v < n; v ++) {
			t2[v] = zint_mod_small_signed(
				gs + v, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t2, t1, p, p0i);
		for (size_t v = 0; v < hn; v ++) {
			yg[v] = mp_montymul(
				mp_montymul(t2[2 * v], t2[2 * v + 1], p, p0i),
				R2, p, p0i);
		}
		yg += hn;
	}
}

/*
 * Compute (f,g) at a specified depth, in RNS+NTT notation.
 * Computed values are stored at the start of the provided tmp[] (slen
 * words per coefficient).
 *
 * This function is for depth < logn_top. For the deepest layer, use
 * make_fg_deepest().
 *
 * RAM USAGE: 3*(2^logn_top)
 */
static void
make_fg_intermediate(const ntru_profile *prof,
	unsigned logn_top,
	const int8_t *restrict f, const int8_t *restrict g,
	unsigned depth, uint32_t *tmp)
{
	make_fg_zero(logn_top, f, g, tmp);
	for (unsigned d = 0; d < depth; d ++) {
		make_fg_step(prof, logn_top, d, tmp);
	}
}

/*
 * Compute (f,g) at the deepest level (i.e. get Res(f,X^n+1) and
 * Res(g,X^n+1)). Intermediate (f,g) values (below the save threshold)
 * are copied at the end of tmp (of size save_off words).
 *
 * If f is not invertible modulo X^n+1 and modulo p = 2147473409, then
 * this function returns 0 (the resultants and intermediate values are
 * still computed). Otherwise it returns 1.
 */
static int
make_fg_deepest(const ntru_profile *prof,
	unsigned logn_top,
	const int8_t *restrict f, const int8_t *restrict g,
	uint32_t *tmp, size_t sav_off)
{
	make_fg_zero(logn_top, f, g, tmp);
	int r = 1;

	/*
	 * f is now in RNS+NTT, so we can test its invertibility (mod p)
	 * by simply checking that all NTT coefficients are non-zero.
	 * (If f is invertible then recover_G() always works.)
	 */
	size_t n = (size_t)1 << logn_top;
	uint32_t b = 0;
	for (size_t u = 0; u < n; u ++) {
		b |= tmp[u] - 1;
	}
	r = (int)(1 - (b >> 31));

	for (unsigned d = 0; d < logn_top; d ++) {
		make_fg_step(prof, logn_top, d, tmp);

		/*
		 * make_fg_step() computes the (f,g) for depth d+1; we
		 * save that value if d+1 is at least at the save
		 * threshold, but is not the deepest level.
		 */
		unsigned d2 = d + 1;
		if (d2 < logn_top && d2 >= prof->min_save_fg[logn_top]) {
			size_t slen = prof->max_bl_small[d2];
			size_t fglen = slen << (logn_top + 1 - d2);
			sav_off -= fglen;
			memmove(tmp + sav_off, tmp, fglen * sizeof *tmp);
		}
	}
	return r;
}

/* Error code: no error (so far) */
#define SOLVE_OK           0

/* Error code: GCD(Res(f,X^n+1), Res(g,X^n+1)) != 1 */
#define SOLVE_ERR_GCD      -1

/* Error code: reduction error (NTRU equation no longer fulfilled) */
#define SOLVE_ERR_REDUCE   -2

/* Error code: output (F,G) coefficients are off-limits */
#define SOLVE_ERR_LIMIT    -3

/*
 * Solve the NTRU equation at the deepest level. This computes the
 * integers F and G such that Res(f,X^n+1)*G - Res(g,X^n+1)*F = q.
 * The two integers are written into tmp[].
 *
 * Returned value: 0 on success, a negative error code otherwise.
 *
 * RAM USAGE: max(3*(2^logn_top), 8*max_bl_small[depth])
 */
static int
solve_NTRU_deepest(const ntru_profile *prof,
	unsigned logn_top, const int8_t *restrict f, const int8_t *restrict g,
	uint32_t *tmp)
{
	/*
	 * Get (f,g) at the deepest level (i.e. Res(f,X^n+1) and Res(g,X^n+1)).
	 * Obtained (f,g) are in RNS+NTT (since degree n = 1, this is
	 * equivalent to RNS).
	 */
	if (!make_fg_deepest(prof, logn_top, f, g,
		tmp, (size_t)6 << logn_top))
	{
		return SOLVE_ERR_GCD;
	}
	/* obsolete
	kg_stats_set_max_size(logn_top, logn_top, 3 * ((size_t)1 << logn_top));
	*/

	/*
	 * Reorganize memory:
	 *    Fp   output F (len)
	 *    Gp   output G (len)
	 *    fp   Res(f,X^n+1) (len)
	 *    gp   Res(g,X^n+1) (len)
	 *    t1   rest of temporary
	 */
	size_t len = prof->max_bl_small[logn_top];
	uint32_t *Fp = tmp;
	uint32_t *Gp = Fp + len;
	uint32_t *fp = Gp + len;
	uint32_t *gp = fp + len;
	uint32_t *t1 = gp + len;
	memmove(fp, tmp, 2 * len * sizeof *tmp);
	/* obsolete
	kg_stats_set_max_size(logn_top, logn_top, 8 * len);
	*/

	/*
	 * Convert back the resultants into plain integers.
	 */
	zint_rebuild_CRT(fp, len, 1, 2, 0, t1);

	/*
	 * Apply the binary GCD to get a solution (F,G) such that
	 * f*G - g*F = 1.
	 */
	if (!zint_bezout(Gp, Fp, fp, gp, len, t1)) {
		return SOLVE_ERR_GCD;
	}

	/*
	 * Multiply the obtained (F,G) by q to get a proper solution
	 * f*G - g*F = q.
	 */
	if (prof->q != 1) {
		if (zint_mul_small(Fp, len, prof->q) != 0
			|| zint_mul_small(Gp, len, prof->q) != 0)
		{
			return SOLVE_ERR_REDUCE;
		}
	}

	return SOLVE_OK;
}

/*
 * We use poly_sub_scaled() when log(n) < MIN_LOGN_FGNTT, and
 * poly_sub_scaled_ntt() when log(n) >= MIN_LOGN_FGNTT. The NTT variant
 * is faster at large degrees, but not at small degrees.
 */
#define MIN_LOGN_FGNTT   4
#if NTRUGEN_AVX2
/*
 * The AVX2 implementation requires MIN_LOGN_FGNTT >= 3
 */
#if MIN_LOGN_FGNTT < 3
#error Incorrect MIN_LOGN_FGNTT value
#endif
#endif // NTRUGEN_AVX2

/*
 * Solving the NTRU equation, intermediate level.
 * Input is (F,G) from one level deeper (half-degree), in plain
 * representation, at the start of tmp[]; output is (F,G) from this
 * level, written at the start of tmp[].
 *
 * Returned value: 0 on success, a negative error code otherwise.
 */
TARGET_AVX2
static int
solve_NTRU_intermediate(const ntru_profile *restrict prof,
	unsigned logn_top,
	const int8_t *restrict f, const int8_t *restrict g,
	unsigned depth, uint32_t *restrict tmp)
{
	/*
	 * MAX SIZE:
	 *    input: 2 * hn * max_bl_small[depth + 1]
	 */

	unsigned logn = logn_top - depth;
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;

	/*
	 * slen   size for (f,g) at this level (also output (F,G))
	 * llen   size for unreduced (F,G) at this level
	 * dlen   size for input (F,G) from deeper level
	 * Note: we always have llen >= dlen (constraint enforced in profiles)
	 */
	size_t slen = prof->max_bl_small[depth];
	size_t llen = prof->max_bl_large[depth];
	size_t dlen = prof->max_bl_small[depth + 1];

	/*
	 * Fd   F from deeper level (dlen*hn)
	 * Gd   G from deeper level (dlen*hn)
	 * ft   f from this level (slen*n)
	 * gt   g from this level (slen*n)
	 */
	uint32_t *Fd = tmp;
	uint32_t *Gd = Fd + dlen * hn;
	uint32_t *fgt = Gd + dlen * hn;

	/*
	 * Get (f,g) for this level (in RNS+NTT).
	 */
	if (depth < prof->min_save_fg[logn_top]) {
		make_fg_intermediate(prof, logn_top, f, g, depth, fgt);
	} else {
		uint32_t *sav_fg = tmp + ((size_t)6 << logn_top);
		for (unsigned d = prof->min_save_fg[logn_top];
			d <= depth; d ++)
		{
			sav_fg -= prof->max_bl_small[d] << (logn_top + 1 - d);
		}
		memmove(fgt, sav_fg, 2 * slen * n * sizeof *fgt);
	}
	/* obsolete
	kg_stats_set_max_size(logn_top, depth,
		2 * dlen * hn + 3 * ((size_t)1 << logn_top));
	*/

	/*
	 * Move buffers so that we have room for the unreduced (F,G) at
	 * this level.
	 *   Ft   F from this level (unreduced) (llen*n)
	 *   Gt   G from this level (unreduced) (llen*n)
	 *   ft   f from this level (slen*n)
	 *   gt   g from this level (slen*n)
	 *   Fd   F from deeper level (dlen*hn)
	 *   Gd   G from deeper level (dlen*hn)
	 */
	uint32_t *Ft = tmp;
	uint32_t *Gt = Ft + llen * n;
	uint32_t *ft = Gt + llen * n;
	uint32_t *gt = ft + slen * n;
	Fd = gt + slen * n;
	Gd = Fd + dlen * hn;
	uint32_t *t1 = Gd + dlen * hn;
	memmove(ft, fgt, 2 * n * slen * sizeof *ft);
	memmove(Fd, tmp, 2 * hn * dlen * sizeof *tmp);
	/* obsolete
	kg_stats_set_max_size(logn_top, depth,
		2 * llen * n + 2 * slen * n + 2 * dlen * hn);
	*/

	/*
	 * Convert Fd and Gd to RNS, with output temporarily stored
	 * in (Ft, Gt). Fd and Gd have degree hn only; we store the
	 * values for each modulus p in the _last_ hn slots of the
	 * n-word line for that modulus.
	 */
	for (size_t u = 0; u < llen; u ++) {
		uint32_t p = PRIMES[u].p;
		uint32_t p0i = PRIMES[u].p0i;
		uint32_t R2 = PRIMES[u].R2;
		uint32_t Rx = mp_Rx31((unsigned)dlen, p, p0i, R2);
		uint32_t *xt = Ft + u * n + hn;
		uint32_t *yt = Gt + u * n + hn;
#if NTRUGEN_AVX2
		if (logn >= 4) {
			__m256i yp = _mm256_set1_epi32(p);
			__m256i yp0i = _mm256_set1_epi32(p0i);
			__m256i yR2 = _mm256_set1_epi32(R2);
			__m256i yRx = _mm256_set1_epi32(Rx);
			for (size_t v = 0; v < hn; v += 8) {
				_mm256_storeu_si256((__m256i *)(xt + v),
					zint_mod_small_signed_x8(Fd + v, dlen,
						hn, yp, yp0i, yR2, yRx));
				_mm256_storeu_si256((__m256i *)(yt + v),
					zint_mod_small_signed_x8(Gd + v, dlen,
						hn, yp, yp0i, yR2, yRx));
			}
		} else {
			for (size_t v = 0; v < hn; v ++) {
				xt[v] = zint_mod_small_signed(Fd + v, dlen, hn,
					p, p0i, R2, Rx);
				yt[v] = zint_mod_small_signed(Gd + v, dlen, hn,
					p, p0i, R2, Rx);
			}
		}
#else // NTRUGEN_AVX2
		for (size_t v = 0; v < hn; v ++) {
			xt[v] = zint_mod_small_signed(Fd + v, dlen, hn,
				p, p0i, R2, Rx);
			yt[v] = zint_mod_small_signed(Gd + v, dlen, hn,
				p, p0i, R2, Rx);
		}
#endif // NTRUGEN_AVX2
	}

	/*
	 * Fd and Gd are no longer needed.
	 */
	t1 = Fd;
	/* obsolete
	kg_stats_set_max_size(logn_top, depth,
		2 * llen * n + 2 * slen * n + llen);
	kg_stats_set_max_size(logn_top, depth,
		2 * llen * n + 2 * slen * n + 4 * n);
	*/

	/*
	 * Compute (F,G) (unreduced) modulo sufficiently many small primes.
	 * We also un-NTT (f,g) as we go; when slen primes have been
	 * processed, we obtain (f,g) in RNS, and we apply the CRT to
	 * get (f,g) in plain representation.
	 */
	for (size_t u = 0; u < llen; u ++) {
		/*
		 * If we have processed exactly slen primes, then (f,g)
		 * are in RNS, and we can rebuild them.
		 */
		if (u == slen) {
			zint_rebuild_CRT(ft, slen, n, 2, 1, t1);
		}

		uint32_t p = PRIMES[u].p;
		uint32_t p0i = PRIMES[u].p0i;
		uint32_t R2 = PRIMES[u].R2;

		/*
		 * Memory layout: we keep Ft, Gt, ft and gt; we append:
		 *   gm    NTT support (n)
		 *   igm   iNTT support (n)
		 *   fx    temporary f mod p (NTT) (n)
		 *   gx    temporary g mod p (NTT) (n)
		 */
		uint32_t *gm = t1;
		uint32_t *igm = gm + n;
		uint32_t *fx = igm + n;
		uint32_t *gx = fx + n;
		mp_mkgmigm(logn, gm, igm, PRIMES[u].g, PRIMES[u].ig, p, p0i);
		if (u < slen) {
			memcpy(fx, ft + u * n, n * sizeof *fx);
			memcpy(gx, gt + u * n, n * sizeof *gx);
			mp_iNTT(logn, ft + u * n, igm, p, p0i);
			mp_iNTT(logn, gt + u * n, igm, p, p0i);
		} else {
			uint32_t Rx = mp_Rx31((unsigned)slen, p, p0i, R2);
			for (size_t v = 0; v < n; v ++) {
				fx[v] = zint_mod_small_signed(ft + v, slen, n,
					p, p0i, R2, Rx);
				gx[v] = zint_mod_small_signed(gt + v, slen, n,
					p, p0i, R2, Rx);
			}
			mp_NTT(logn, fx, gm, p, p0i);
			mp_NTT(logn, gx, gm, p, p0i);
		}

		/*
		 * We have (F,G) from deeper level in Ft and Gt, in
		 * RNS. We apply the NTT modulo p.
		 */
		uint32_t *Fe = Ft + u * n;
		uint32_t *Ge = Gt + u * n;
		mp_NTT(logn - 1, Fe + hn, gm, p, p0i);
		mp_NTT(logn - 1, Ge + hn, gm, p, p0i);

		/*
		 * Compute F and G (unreduced) modulo p.
		 */
#if NTRUGEN_AVX2
		if (hn >= 4) {
			__m256i yp = _mm256_set1_epi32(p);
			__m256i yp0i = _mm256_set1_epi32(p0i);
			__m256i yR2 = _mm256_set1_epi32(R2);
			for (size_t v = 0; v < hn; v += 4) {
				__m256i yfa = _mm256_loadu_si256(
					(__m256i *)(fx + (v << 1)));
				__m256i yga = _mm256_loadu_si256(
					(__m256i *)(gx + (v << 1)));
				__m256i yfb = _mm256_srli_epi64(yfa, 32);
				__m256i ygb = _mm256_srli_epi64(yga, 32);
				__m128i xFe = _mm_loadu_si128(
					(__m128i *)(Fe + v + hn));
				__m128i xGe = _mm_loadu_si128(
					(__m128i *)(Ge + v + hn));
				__m256i yFp = _mm256_permute4x64_epi64(
					_mm256_castsi128_si256(xFe), 0x50);
				__m256i yGp = _mm256_permute4x64_epi64(
					_mm256_castsi128_si256(xGe), 0x50);
				yFp = _mm256_shuffle_epi32(yFp, 0x30);
				yGp = _mm256_shuffle_epi32(yGp, 0x30);
				yFp = mp_montymul_x4(yFp, yR2, yp, yp0i);
				yGp = mp_montymul_x4(yGp, yR2, yp, yp0i);
				__m256i yFe0 = mp_montymul_x4(
					ygb, yFp, yp, yp0i);
				__m256i yFe1 = mp_montymul_x4(
					yga, yFp, yp, yp0i);
				__m256i yGe0 = mp_montymul_x4(
					yfb, yGp, yp, yp0i);
				__m256i yGe1 = mp_montymul_x4(
					yfa, yGp, yp, yp0i);
				_mm256_storeu_si256((__m256i *)(Fe + (v << 1)),
					_mm256_or_si256(yFe0,
						_mm256_slli_epi64(yFe1, 32)));
				_mm256_storeu_si256((__m256i *)(Ge + (v << 1)),
					_mm256_or_si256(yGe0,
						_mm256_slli_epi64(yGe1, 32)));
			}
		} else {
			for (size_t v = 0; v < hn; v ++) {
				uint32_t fa = fx[(v << 1) + 0];
				uint32_t fb = fx[(v << 1) + 1];
				uint32_t ga = gx[(v << 1) + 0];
				uint32_t gb = gx[(v << 1) + 1];
				uint32_t mFp = mp_montymul(
					Fe[v + hn], R2, p, p0i);
				uint32_t mGp = mp_montymul(
					Ge[v + hn], R2, p, p0i);
				Fe[(v << 1) + 0] = mp_montymul(gb, mFp, p, p0i);
				Fe[(v << 1) + 1] = mp_montymul(ga, mFp, p, p0i);
				Ge[(v << 1) + 0] = mp_montymul(fb, mGp, p, p0i);
				Ge[(v << 1) + 1] = mp_montymul(fa, mGp, p, p0i);
			}
		}
#else // NTRUGEN_AVX2
		for (size_t v = 0; v < hn; v ++) {
			uint32_t fa = fx[(v << 1) + 0];
			uint32_t fb = fx[(v << 1) + 1];
			uint32_t ga = gx[(v << 1) + 0];
			uint32_t gb = gx[(v << 1) + 1];
			uint32_t mFp = mp_montymul(Fe[v + hn], R2, p, p0i);
			uint32_t mGp = mp_montymul(Ge[v + hn], R2, p, p0i);
			Fe[(v << 1) + 0] = mp_montymul(gb, mFp, p, p0i);
			Fe[(v << 1) + 1] = mp_montymul(ga, mFp, p, p0i);
			Ge[(v << 1) + 0] = mp_montymul(fb, mGp, p, p0i);
			Ge[(v << 1) + 1] = mp_montymul(fa, mGp, p, p0i);
		}
#endif // NTRUGEN_AVX2

		/*
		 * We want the new (F,G) in RNS only (no NTT).
		 */
		mp_iNTT(logn, Fe, igm, p, p0i);
		mp_iNTT(logn, Ge, igm, p, p0i);
	}

	/*
	 * Edge case: if slen == llen, then we have not rebuilt (f,g)
	 * into plain representation yet, so we do it now.
	 */
	if (slen == llen) {
		zint_rebuild_CRT(ft, slen, n, 2, 1, t1);
	}

	/*
	 * We now have the unreduced (F,G) in RNS. We rebuild their
	 * plain representation.
	 */
	zint_rebuild_CRT(Ft, llen, n, 2, 1, t1);

	/*
	 * We now reduce these (F,G) with Babai's nearest plane
	 * algorithm. The reduction conceptually goes as follows:
	 *   k <- round((F*adj(f) + G*adj(g))/(f*adj(f) + g*adj(g)))
	 *   (F, G) <- (F - k*f, G - k*g)
	 * We use fixed-point approximations of (f,g) and (F, G) to get
	 * a value k as a small polynomial with scaling; we then apply
	 * k on the full-width polynomial. Each iteration "shaves" a
	 * a few bits off F and G.
	 *
	 * We apply the process sufficiently many times to reduce (F, G)
	 * to the size of (f, g) with a reasonable probability of success.
	 * Since we want full constant-time processing, the number of
	 * iterations and the accessed slots work on some assumptions on
	 * the sizes of values (sizes have been measured over many samples,
	 * and a margin of 5 times the standard deviation).
	 */

	/*
	 * If depth is at least 2, and we will use the NTT to subtract
	 * k*(f,g) from (F,G), then we will need to convert (f,g) to NTT over
	 * slen+1 words, which requires an extra word to ft and gt.
	 */
	int use_sub_ntt = (depth > 1 && logn >= MIN_LOGN_FGNTT);
	if (use_sub_ntt) {
		memmove(gt + n, gt, n * slen * sizeof *gt);
		gt += n;
		t1 += 2 * n;
	}

	/*
	 * New layout:
	 *   Ft    F from this level (unreduced) (llen*n)
	 *   Gt    G from this level (unreduced) (llen*n)
	 *   ft    f from this level (slen*n) (+n if use_sub_ntt)
	 *   gt    g from this level (slen*n) (+n if use_sub_ntt)
	 *   rt3   (n fxr = 2*n)
	 *   rt4   (n fxr = 2*n)
	 *   rt1   (hn fxr = n)
	 */

	fxr *rt3 = (fxr *)t1;
	fxr *rt4 = rt3 + n;
	fxr *rt1 = rt4 + n;
	/* obsolete
	kg_stats_set_max_size(logn_top, depth,
		2 * llen * n + 2 * slen * n + 5 * n);
	*/

	/*
	 * We consider only the top rlen words of (f,g).
	 */
	size_t rlen = prof->word_win[depth];
	if (rlen > slen) {
		rlen = slen;
	}
	size_t blen = slen - rlen;
	uint32_t *ftb = ft + blen * n;
	uint32_t *gtb = gt + blen * n;
	uint32_t scale_fg = 31 * (uint32_t)blen;
	uint32_t scale_FG = 31 * (uint32_t)llen;

	/*
	 * Convert f and g into fixed-point approximations, in rt3 and rt4,
	 * respectively. They are scaled down by 2^(scale_fg + scale_x).
	 * scale_fg is public (it depends only on the recursion depth), but
	 * scale_x comes from a measurement on the actual values of (f,g) and
	 * is thus secret.
	 *
	 * The value scale_x is adjusted so that the largest coefficient is
	 * close to, but lower than, some limit t (in absolute value). The
	 * limit t is chosen so that f*adj(f) + g*adj(g) does not overflow,
	 * i.e. all coefficients must remain below 2^31.
	 *
	 * Let n be the degree; we know that n <= 2^10. The squared norm
	 * of a polynomial is the sum of the squared norms of the
	 * coefficients, with the squared norm of a complex number being
	 * the product of that number with its complex conjugate. If all
	 * coefficients of f are less than t (in absolute value), then
	 * the squared norm of f is less than n*t^2. The squared norm of
	 * FFT(f) (f in FFT representation) is exactly n times the
	 * squared norm of f, so this leads to n^2*t^2 as a maximum
	 * bound. adj(f) has the same norm as f. This implies that each
	 * complex coefficient of FFT(f) has a maximum squared norm of
	 * n^2*t^2 (with a maximally imbalanced polynomial with all
	 * coefficient but one being zero). The computation of f*adj(f)
	 * exactly is, in FFT representation, the product of each
	 * coefficient with its conjugate; thus, the coefficients of
	 * f*adj(f), in FFT representation, are at most n^2*t^2.
	 *
	 * Since we want the coefficients of f*adj(f)+g*adj(g) not to exceed
	 * 2^31, we need n^2*t^2 <= 2^30, i.e. n*t <= 2^15. We can adjust t
	 * accordingly (called scale_t in the code below). We also need to
	 * take care that t must not exceed scale_x. Approximation of f and
	 * g are extracted with scale scale_fg + scale_x - scale_t, and
	 * later fixed by dividing them by 2^scale_t.
	 */
	uint32_t scale_xf = poly_max_bitlength(logn, ftb, rlen);
	uint32_t scale_xg = poly_max_bitlength(logn, gtb, rlen);
	uint32_t scale_x = scale_xf;
	scale_x ^= (scale_xf ^ scale_xg) & tbmask(scale_xf - scale_xg);
	uint32_t scale_t = 15 - logn;
	scale_t ^= (scale_t ^ scale_x) & tbmask(scale_x - scale_t);
	uint32_t scdiff = scale_x - scale_t;

	poly_big_to_fixed(logn, rt3, ftb, rlen, scdiff);
	poly_big_to_fixed(logn, rt4, gtb, rlen, scdiff);

	/*
	 * Compute adj(f)/(f*adj(f) + g*adj(g)) into rt3 (FFT).
	 * Compute adj(g)/(f*adj(f) + g*adj(g)) into rt4 (FFT).
	 */
	vect_FFT(logn, rt3);
	vect_FFT(logn, rt4);
	vect_norm_fft(logn, rt1, rt3, rt4);
	vect_mul2e(logn, rt3, scale_t);
	vect_mul2e(logn, rt4, scale_t);
	for (size_t u = 0; u < hn; u ++) {
#if NTRUGEN_AVX2
		fxr ni3 = fxr_neg(rt3[u + hn]);
		fxr ni4 = fxr_neg(rt4[u + hn]);
		fxr_div_x4_1(&rt3[u], &ni3, &rt4[u], &ni4, rt1[u]);
		rt3[u + hn] = ni3;
		rt4[u + hn] = ni4;
#else // NTRUGEN_AVX2
		rt3[u] = fxr_div(rt3[u], rt1[u]);
		rt3[u + hn] = fxr_div(fxr_neg(rt3[u + hn]), rt1[u]);
		rt4[u] = fxr_div(rt4[u], rt1[u]);
		rt4[u + hn] = fxr_div(fxr_neg(rt4[u + hn]), rt1[u]);
#endif // NTRUGEN_AVX2
	}

	/*
	 * New layout:
	 *   Ft    F from this level (unreduced) (llen*n)
	 *   Gt    G from this level (unreduced) (llen*n)
	 *   ft    f from this level (slen*n) (+n if use_sub_ntt)
	 *   gt    g from this level (slen*n) (+n if use_sub_ntt)
	 *   rt3   (n fxr = 2*n)
	 *   rt4   (n fxr = 2*n)
	 *   rt1   (n fxr = 2*n)     |   k    (n)
	 *   rt2   (n fxr = 2*n)     |   t2   (3*n)
	 * Exception: at depth == 1, we omit ft and gt:
	 *   Ft    F from this level (unreduced) (llen*n)
	 *   Gt    G from this level (unreduced) (llen*n)
	 *   rt3   (n fxr = 2*n)
	 *   rt4   (n fxr = 2*n)
	 *   rt1   (n fxr = 2*n)     |   k    (n)
	 *   rt2   (n fxr = 2*n)     |   t2   (3*n)
	 */
	if (depth == 1) {
		t1 = ft;
		fxr *nrt3 = (fxr *)t1;
		memmove(nrt3, rt3, 2 * n * sizeof *rt3);
		rt3 = nrt3;
		rt4 = rt3 + n;
		rt1 = rt4 + n;
	}
	int32_t *k = (int32_t *)rt1;
	uint32_t *t2 = (uint32_t *)(k + n);
	fxr *rt2 = (fxr *)t2;
	if (rt2 < (rt1 + n)) {
		rt2 = rt1 + n;
	}
	/* obsolete
	if (depth == 1) {
		kg_stats_set_max_size(logn_top, depth,
			2 * llen * n + 8 * n);
	} else {
		kg_stats_set_max_size(logn_top, depth,
			2 * llen * n + 2 * slen * n + 8 * n);
		kg_stats_set_max_size(logn_top, depth,
			2 * llen * n + 2 * slen * n + 4 * n + n
			+ 2 * n + (slen + 1) * n + n);
	}
	*/

	/*
	 * If we are going to use poly_sub_scaled_ntt(), then we convert
	 * f and g to the NTT representation. Since poly_sub_scaled_ntt()
	 * itself will use more than n*(slen+2) words in t2[], we can do
	 * the same here.
	 */
	if (use_sub_ntt) {
		uint32_t *gm = t2;
		uint32_t *tn = gm + n;
		for (size_t u = 0; u <= slen; u ++) {
			uint32_t p = PRIMES[u].p;
			uint32_t p0i = PRIMES[u].p0i;
			uint32_t R2 = PRIMES[u].R2;
			uint32_t Rx = mp_Rx31((unsigned)slen, p, p0i, R2);
			mp_mkgm(logn, gm, PRIMES[u].g, p, p0i);
			for (size_t v = 0; v < n; v ++) {
				tn[v] = zint_mod_small_signed(
					ft + v, slen, n, p, p0i, R2, Rx);
			}
			mp_NTT(logn, tn, gm, p, p0i);
			tn += n;
		}
		tn = gm + n;
		memmove(ft, tn, (slen + 1) * n * sizeof *tn);
		for (size_t u = 0; u <= slen; u ++) {
			uint32_t p = PRIMES[u].p;
			uint32_t p0i = PRIMES[u].p0i;
			uint32_t R2 = PRIMES[u].R2;
			uint32_t Rx = mp_Rx31((unsigned)slen, p, p0i, R2);
			mp_mkgm(logn, gm, PRIMES[u].g, p, p0i);
			for (size_t v = 0; v < n; v ++) {
				tn[v] = zint_mod_small_signed(
					gt + v, slen, n, p, p0i, R2, Rx);
			}
			mp_NTT(logn, tn, gm, p, p0i);
			tn += n;
		}
		tn = gm + n;
		memmove(gt, tn, (slen + 1) * n * sizeof *tn);
	}

	/*
	 * Reduce F and G repeatedly.
	 */
	size_t FGlen = llen;
	for (;;) {
		/*
		 * Convert the current F and G into fixed-point. We want
		 * to apply scaling scale_FG + scale_x.
		 */
		uint32_t tlen, toff;
		DIVREM31(tlen, toff, scale_FG);
		poly_big_to_fixed(logn, rt1,
			Ft + tlen * n, FGlen - tlen, scale_x + toff);
		poly_big_to_fixed(logn, rt2,
			Gt + tlen * n, FGlen - tlen, scale_x + toff);

		/*
		 * rt2 <- (F*adj(f) + G*adj(g)) / (f*adj(f) + g*adj(g))
		 */
		vect_FFT(logn, rt1);
		vect_FFT(logn, rt2);
		vect_mul_fft(logn, rt1, rt3);
		vect_mul_fft(logn, rt2, rt4);
		vect_add(logn, rt2, rt1);
		vect_iFFT(logn, rt2);

		/*
		 * k <- round(rt2)
		 */
		for (size_t u = 0; u < n; u ++) {
			k[u] = fxr_round(rt2[u]);
		}

		/*
		 * (f,g) are scaled by scale_fg + scale_x
		 * (F,G) are scaled by scale_FG + scale_x
		 * Thus, k is scaled by scale_FG - scale_fg, which is public.
		 */
		uint32_t scale_k = scale_FG - scale_fg;
		if (depth == 1) {
			poly_sub_kfg_scaled_depth1(logn_top, Ft, Gt, FGlen,
				(uint32_t *)k, scale_k, f, g, t2);
		} else if (use_sub_ntt) {
			poly_sub_scaled_ntt(logn, Ft, FGlen, ft, slen,
				k, scale_k, t2);
			poly_sub_scaled_ntt(logn, Gt, FGlen, gt, slen,
				k, scale_k, t2);
		} else {
			poly_sub_scaled(logn, Ft, FGlen, ft, slen, k, scale_k);
			poly_sub_scaled(logn, Gt, FGlen, gt, slen, k, scale_k);
		}

		/*
		 * We now assume that F and G have shrunk by at least
		 * reduce_bits (profile-dependent). We adjust FGlen accordinly.
		 */
		if (scale_FG <= scale_fg) {
			break;
		}
		if (scale_FG <= (scale_fg + prof->reduce_bits)) {
			scale_FG = scale_fg;
		} else {
			scale_FG -= prof->reduce_bits;
		}
		while (FGlen > slen
			&& 31 * (FGlen - slen) > scale_FG - scale_fg + 30)
		{
			FGlen --;
		}
	}

	/*
	 * Output F is already in the right place; G is in Gt, and must be
	 * moved back a bit.
	 */
	memmove(tmp + slen * n, Gt, slen * n * sizeof *tmp);
	Gt = tmp + slen * n;

	/*
	 * Reduction is done. We test the current solution modulo a single
	 * prime.
	 * Exception: we cannot do that if depth == 1, since in that case
	 * we did not keep (ft,gt). Reduction errors rarely occur at this
	 * stage, so we can omit that test (depth-0 test will cover it).
	 *
	 * If use_sub_ntt != 0, then ft and gt are already in NTT
	 * representation.
	 */
	if (depth == 1) {
		return SOLVE_OK;
	}

	t2 = t1 + n;
	uint32_t *t3 = t2 + n;
	uint32_t *t4 = t3 + n;
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	uint32_t R2 = PRIMES[0].R2;
	uint32_t Rx = mp_Rx31(slen, p, p0i, R2);
	mp_mkgm(logn, t4, PRIMES[0].g, p, p0i);
	if (use_sub_ntt) {
		t1 = ft;
		for (size_t u = 0; u < n; u ++) {
			t2[u] = zint_mod_small_signed(
				Gt + u, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t2, t4, p, p0i);
	} else {
		for (size_t u = 0; u < n; u ++) {
			t1[u] = zint_mod_small_signed(
				ft + u, slen, n, p, p0i, R2, Rx);
			t2[u] = zint_mod_small_signed(
				Gt + u, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t1, t4, p, p0i);
		mp_NTT(logn, t2, t4, p, p0i);
	}
#if NTRUGEN_AVX2
	if (n >= 8) {
		__m256i yp = _mm256_set1_epi32(p);
		__m256i yp0i = _mm256_set1_epi32(p0i);
		for (size_t u = 0; u < n; u += 8) {
			__m256i y1 = _mm256_loadu_si256((__m256i *)(t1 + u));
			__m256i y2 = _mm256_loadu_si256((__m256i *)(t2 + u));
			__m256i y3 = mp_montymul_x8(y1, y2, yp, yp0i);
			_mm256_storeu_si256((__m256i *)(t3 + u), y3);
		}
	} else {
		for (size_t u = 0; u < n; u ++) {
			t3[u] = mp_montymul(t1[u], t2[u], p, p0i);
		}
	}
#else // NTRUGEN_AVX2
	for (size_t u = 0; u < n; u ++) {
		t3[u] = mp_montymul(t1[u], t2[u], p, p0i);
	}
#endif // NTRUGEN_AVX2
	if (use_sub_ntt) {
		t1 = gt;
		for (size_t u = 0; u < n; u ++) {
			t2[u] = zint_mod_small_signed(
				Ft + u, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t2, t4, p, p0i);
	} else {
		for (size_t u = 0; u < n; u ++) {
			t1[u] = zint_mod_small_signed(
				gt + u, slen, n, p, p0i, R2, Rx);
			t2[u] = zint_mod_small_signed(
				Ft + u, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t1, t4, p, p0i);
		mp_NTT(logn, t2, t4, p, p0i);
	}
	uint32_t rv = mp_montymul(prof->q, 1, p, p0i);
#if NTRUGEN_AVX2
	if (n >= 8) {
		__m256i yp = _mm256_set1_epi32(p);
		__m256i yp0i = _mm256_set1_epi32(p0i);
		__m256i yrv = _mm256_set1_epi32(rv);
		for (size_t u = 0; u < n; u += 8) {
			__m256i y1 = _mm256_loadu_si256((__m256i *)(t1 + u));
			__m256i y2 = _mm256_loadu_si256((__m256i *)(t2 + u));
			__m256i y3 = _mm256_loadu_si256((__m256i *)(t3 + u));
			__m256i yx = mp_sub_x8(y3,
				mp_montymul_x8(y1, y2, yp, yp0i), yp);
			if ((uint32_t)_mm256_movemask_epi8(
				_mm256_cmpeq_epi32(yx, yrv)) != 0xFFFFFFFF)
			{
				return SOLVE_ERR_REDUCE;
			}
		}
	} else {
		for (size_t u = 0; u < n; u ++) {
			uint32_t x = mp_montymul(t1[u], t2[u], p, p0i);
			if (mp_sub(t3[u], x, p) != rv) {
				return SOLVE_ERR_REDUCE;
			}
		}
	}
#else // NTRUGEN_AVX2
	for (size_t u = 0; u < n; u ++) {
		uint32_t x = mp_montymul(t1[u], t2[u], p, p0i);
		if (mp_sub(t3[u], x, p) != rv) {
			return SOLVE_ERR_REDUCE;
		}
	}
#endif // NTRUGEN_AVX2

	return SOLVE_OK;
}

/*
 * Solving the NTRU equation, top recursion level. This is a specialized
 * variant for solve_NTRU_intermediate() with depth == 0, for lower RAM
 * usage and faster operation.
 *
 * Returned value: 0 on success, a negative error code otherwise.
 */
TARGET_AVX2
static int
solve_NTRU_depth0(const ntru_profile *restrict prof,
	unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g,
	uint32_t *restrict tmp)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;

	/*
	 * At depth 0, all values fit on 30 bits, so we work with a
	 * single modulus p.
	 */
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	uint32_t R2 = PRIMES[0].R2;

	/*
	 * Buffer layout:
	 *   Fd   F from upper level (hn)
	 *   Gd   G from upper level (hn)
	 *   ft   f (n)
	 *   gt   g (n)
	 *   gm   helper for NTT
	 */
	uint32_t *Fd = tmp;
	uint32_t *Gd = Fd + hn;
	uint32_t *ft = Gd + hn;
	uint32_t *gt = ft + n;
	uint32_t *gm = gt + n;
	/* obsolete
	kg_stats_set_max_size(logn, 0, 2 * hn + 3 * n);
	*/

	/*
	 * Load f and g, and convert to RNS+NTT.
	 */
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	poly_mp_set_small(logn, ft, f, p);
	poly_mp_set_small(logn, gt, g, p);
	mp_NTT(logn, ft, gm, p, p0i);
	mp_NTT(logn, gt, gm, p, p0i);

	/*
	 * Convert Fd and Gd to RNS+NTT.
	 */
	poly_mp_set(logn - 1, Fd, p);
	poly_mp_set(logn - 1, Gd, p);
	mp_NTT(logn - 1, Fd, gm, p, p0i);
	mp_NTT(logn - 1, Gd, gm, p, p0i);

	/*
	 * Build the unreduced (F,G) into ft and gt.
	 */
#if NTRUGEN_AVX2
	if (hn >= 4) {
		__m256i yp = _mm256_set1_epi32(p);
		__m256i yp0i = _mm256_set1_epi32(p0i);
		__m256i yR2 = _mm256_set1_epi32(R2);
		for (size_t v = 0; v < hn; v += 4) {
			__m256i yfa = _mm256_loadu_si256(
				(__m256i *)(ft + (v << 1)));
			__m256i yga = _mm256_loadu_si256(
				(__m256i *)(gt + (v << 1)));
			__m256i yfb = _mm256_srli_epi64(yfa, 32);
			__m256i ygb = _mm256_srli_epi64(yga, 32);
			__m128i xFd = _mm_loadu_si128((__m128i *)(Fd + v));
			__m128i xGd = _mm_loadu_si128((__m128i *)(Gd + v));
			__m256i yFd = _mm256_permute4x64_epi64(
				_mm256_castsi128_si256(xFd), 0x50);
			__m256i yGd = _mm256_permute4x64_epi64(
				_mm256_castsi128_si256(xGd), 0x50);
			yFd = _mm256_shuffle_epi32(yFd, 0x30);
			yGd = _mm256_shuffle_epi32(yGd, 0x30);
			yFd = mp_montymul_x4(yFd, yR2, yp, yp0i);
			yGd = mp_montymul_x4(yGd, yR2, yp, yp0i);
			__m256i yFe0 = mp_montymul_x4(ygb, yFd, yp, yp0i);
			__m256i yFe1 = mp_montymul_x4(yga, yFd, yp, yp0i);
			__m256i yGe0 = mp_montymul_x4(yfb, yGd, yp, yp0i);
			__m256i yGe1 = mp_montymul_x4(yfa, yGd, yp, yp0i);
			_mm256_storeu_si256((__m256i *)(ft + (v << 1)),
				_mm256_or_si256(yFe0,
					_mm256_slli_epi64(yFe1, 32)));
			_mm256_storeu_si256((__m256i *)(gt + (v << 1)),
				_mm256_or_si256(yGe0,
					_mm256_slli_epi64(yGe1, 32)));
		}
	} else {
		for (size_t v = 0; v < hn; v ++) {
			uint32_t fa = ft[(v << 1) + 0];
			uint32_t fb = ft[(v << 1) + 1];
			uint32_t ga = gt[(v << 1) + 0];
			uint32_t gb = gt[(v << 1) + 1];
			uint32_t mFd = mp_montymul(Fd[v], R2, p, p0i);
			uint32_t mGd = mp_montymul(Gd[v], R2, p, p0i);
			ft[(v << 1) + 0] = mp_montymul(gb, mFd, p, p0i);
			ft[(v << 1) + 1] = mp_montymul(ga, mFd, p, p0i);
			gt[(v << 1) + 0] = mp_montymul(fb, mGd, p, p0i);
			gt[(v << 1) + 1] = mp_montymul(fa, mGd, p, p0i);
		}
	}
#else // NTRUGEN_AVX2
	for (size_t u = 0; u < hn; u ++) {
		uint32_t fa = ft[(u << 1) + 0];
		uint32_t fb = ft[(u << 1) + 1];
		uint32_t ga = gt[(u << 1) + 0];
		uint32_t gb = gt[(u << 1) + 1];
		uint32_t mFd = mp_montymul(Fd[u], R2, p, p0i);
		uint32_t mGd = mp_montymul(Gd[u], R2, p, p0i);
		ft[(u << 1) + 0] = mp_montymul(gb, mFd, p, p0i);
		ft[(u << 1) + 1] = mp_montymul(ga, mFd, p, p0i);
		gt[(u << 1) + 0] = mp_montymul(fb, mGd, p, p0i);
		gt[(u << 1) + 1] = mp_montymul(fa, mGd, p, p0i);
	}
#endif // NTRUGEN_AVX2

	/*
	 * Reorganize buffers:
	 *   Fp   unreduced F (RNS+NTT) (n)
	 *   Gp   unreduced G (RNS+NTT) (n)
	 *   t1   free (n)
	 *   t2   NTT support (gm) (n)
	 *   t3   free (n)
	 *   t4   free (n)
	 */
	uint32_t *Fp = tmp;
	uint32_t *Gp = Fp + n;
	uint32_t *t1 = Gp + n;
	uint32_t *t2 = t1 + n;  /* alias on gm */
	uint32_t *t3 = t2 + n;
	uint32_t *t4 = t3 + n;
	memmove(Fp, ft, 2 * n * sizeof *ft);
	/* obsolete
	kg_stats_set_max_size(logn, 0, 6 * n);
	*/

	/*
	 * Working modulo p (using the NTT), we compute:
	 *    t1 <- F*adj(f) + G*adj(g)
	 *    t2 <- f*adj(f) + g*adj(g)
	 */

	/*
	 * t4 <- f (RNS+NTT)
	 */
	poly_mp_set_small(logn, t4, f, p);
	mp_NTT(logn, t4, gm, p, p0i);

	/*
	 * t1 <- F*adj(f) (RNS+NTT)
	 * t3 <- f*adj(f) (RNS+NTT)
	 */
	for (size_t u = 0; u < n; u ++) {
		uint32_t w = mp_montymul(t4[(n - 1) - u], R2, p, p0i);
		t1[u] = mp_montymul(w, Fp[u], p, p0i);
		t3[u] = mp_montymul(w, t4[u], p, p0i);
	}

	/*
	 * t4 <- g (RNS+NTT)
	 */
	poly_mp_set_small(logn, t4, g, p);
	mp_NTT(logn, t4, gm, p, p0i);

	/*
	 * t1 <- t1 + G*adj(g)
	 * t3 <- t3 + g*adj(g)
	 */
	for (size_t u = 0; u < n; u ++) {
		uint32_t w = mp_montymul(t4[(n - 1) - u], R2, p, p0i);
		t1[u] = mp_add(t1[u], mp_montymul(w, Gp[u], p, p0i), p);
		t3[u] = mp_add(t3[u], mp_montymul(w, t4[u], p, p0i), p);
	}

	/*
	 * Convert back F*adj(f) + G*adj(g) and f*adj(f) + g*adj(g) to
	 * plain representation, and move f*adj(f) + g*adj(g) to t2.
	 */
	mp_mkigm(logn, t4, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t1, t4, p, p0i);
	mp_iNTT(logn, t3, t4, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		/*
		 * NOTE: no truncature to 31 bits.
		 */
		t1[u] = (uint32_t)mp_norm(t1[u], p);
		t2[u] = (uint32_t)mp_norm(t3[u], p);
	}

	/*
	 * Buffer contents:
	 *   Fp   unreduced F (RNS+NTT) (n)
	 *   Gp   unreduced G (RNS+NTT) (n)
	 *   t1   F*adj(f) + G*adj(g) (plain, 32-bit) (n)
	 *   t2   f*adj(f) + g*adj(g) (plain, 32-bit) (n)
	 */

	/*
	 * We need to divide t1 by t2, and round the result. We convert
	 * them to FFT representation, downscaled by 2^10 (to avoid overflows).
	 * We first convert f*adj(f) + g*adj(g), which is auto-adjoint;
	 * thus, its FFT representation only has half-size.
	 */
	fxr *rt3 = (fxr *)t3;
	for (size_t u = 0; u < n; u ++) {
		uint64_t x = (uint64_t)*(int32_t *)&t2[u] << 22;
		rt3[u] = fxr_of_scaled32(x);
	}
	vect_FFT(logn, rt3);
	fxr *rt2 = (fxr *)t2;
	memmove(rt2, rt3, hn * sizeof *rt3);
	rt3 = rt2 + hn;
	/* obsolete
	kg_stats_set_max_size(logn, 0, 6 * n);
	*/

	/*
	 * Buffer contents:
	 *   Fp    unreduced F (RNS+NTT) (n)
	 *   Gp    unreduced G (RNS+NTT) (n)
	 *   t1    F*adj(f) + G*adj(g) (plain, 32-bit) (n)
	 *   rt2   f*adj(f) + g*adj(g) (FFT, auto-ajdoint) (hn fxr values = n)
	 *   rt3   free (n fxr values = 2*n)
	 */

	/*
	 * Convert F*adj(f) + G*adj(g) to FFT (scaled by 2^10) (into rt3).
	 */
	for (size_t u = 0; u < n; u ++) {
		uint64_t x = (uint64_t)*(int32_t *)&t1[u] << 22;
		rt3[u] = fxr_of_scaled32(x);
	}
	vect_FFT(logn, rt3);

	/*
	 * Divide F*adj(f) + G*adj(g) by f*adj(f) + g*adj(g) and round
	 * the result into t1, with conversion to RNS.
	 */
	vect_div_autoadj_fft(logn, rt3, rt2);
	vect_iFFT(logn, rt3);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mp_set(fxr_round(rt3[u]), p);
	}

	/*
	 * Buffer contents:
	 *   Fp    unreduced F (RNS+NTT) (n)
	 *   Gp    unreduced G (RNS+NTT) (n)
	 *   t1    k (RNS) (n)
	 *   t2    free (n)
	 *   t3    free (n)
	 *   t4    free (n)
	 */

	/*
	 * Convert k to RNS+NTT+Montgomery.
	 */
	mp_mkgm(logn, t4, PRIMES[0].g, p, p0i);
	mp_NTT(logn, t1, t4, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mp_montymul(t1[u], R2, p, p0i);
	}

	/*
	 * Subtract k*f from F and k*g from G.
	 * We also compute f*G - g*F (in RNS+NTT) to check that the solution
	 * is correct.
	 */
	for (size_t u = 0; u < n; u ++) {
		t2[u] = mp_set(f[u], p);
		t3[u] = mp_set(g[u], p);
	}
	mp_NTT(logn, t2, t4, p, p0i);
	mp_NTT(logn, t3, t4, p, p0i);
	uint32_t rv = mp_montymul(prof->q, 1, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		Fp[u] = mp_sub(Fp[u], mp_montymul(t1[u], t2[u], p, p0i), p);
		Gp[u] = mp_sub(Gp[u], mp_montymul(t1[u], t3[u], p, p0i), p);
		uint32_t x = mp_sub(
			mp_montymul(t2[u], Gp[u], p, p0i),
			mp_montymul(t3[u], Fp[u], p, p0i), p);
		if (x != rv) {
			return SOLVE_ERR_REDUCE;
		}
	}

	/*
	 * Convert back F and G into normal representation.
	 */
	mp_mkigm(logn, t4, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, Fp, t4, p, p0i);
	mp_iNTT(logn, Gp, t4, p, p0i);
	poly_mp_norm(logn, Fp, p);
	poly_mp_norm(logn, Gp, p);

	return SOLVE_OK;
}

#if 0
/* unused */
/*
 * Verify that the given f, g, F and G fulfill the NTRU equation.
 * Returned value is 1 on success, 0 on error.
 *
 * RAM USAGE: 4*n words
 */
static int
verify_NTRU(const ntru_profile *restrict prof, unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	uint32_t *tmp)
{
	size_t n = (size_t)1 << logn;
	uint32_t *t1 = tmp;
	uint32_t *t2 = t1 + n;
	uint32_t *t3 = t2 + n;
	uint32_t *t4 = t3 + n;
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	mp_mkgm(logn, t4, PRIMES[0].g, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mp_set(f[u], p);
		t2[u] = mp_set(G[u], p);
	}
	mp_NTT(logn, t1, t4, p, p0i);
	mp_NTT(logn, t2, t4, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		t3[u] = mp_montymul(t1[u], t2[u], p, p0i);
	}
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mp_set(g[u], p);
		t2[u] = mp_set(F[u], p);
	}
	mp_NTT(logn, t1, t4, p, p0i);
	mp_NTT(logn, t2, t4, p, p0i);
	uint32_t rv = mp_montymul(prof->q, 1, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		uint32_t x = mp_montymul(t1[u], t2[u], p, p0i);
		if (mp_sub(t3[u], x, p) != rv) {
			return 0;
		}
	}
	return 1;
}
#endif

/* see ng_inner.h */
int
solve_NTRU(const ntru_profile *restrict prof, unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g, uint32_t *tmp)
{
	size_t n = (size_t)1 << logn;

	int err = solve_NTRU_deepest(prof, logn, f, g, tmp);
	if (err != SOLVE_OK) {
		return err;
	}
	unsigned depth = logn;
	while (depth -- > 1) {
		err = solve_NTRU_intermediate(prof, logn, f, g, depth, tmp);
		if (err != SOLVE_OK) {
			return err;
		}
	}
	err = solve_NTRU_depth0(prof, logn, f, g, tmp);
	if (err != SOLVE_OK) {
		return err;
	}

	/*
	 * F and G are at the start of tmp[] (plain, 31 bits per value).
	 * We need to convert them to 8-bit representation, and check
	 * that they are within the expected range.
	 */
	int8_t *F = (int8_t *)(tmp + 2 * n);
	int8_t *G = F + n;
	int lim = prof->coeff_FG_limit[logn];
	if (!poly_big_to_small(logn, F, tmp, lim)) {
		return SOLVE_ERR_LIMIT;
	}
	if (!poly_big_to_small(logn, G, tmp + n, lim)) {
		return SOLVE_ERR_LIMIT;
	}
	memmove(tmp, F, 2 * n);

	return SOLVE_OK;
}

/* see ng_inner.h */
int
recover_G(unsigned logn, int32_t q, uint32_t ulim,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, uint32_t *restrict tmp)
{
	size_t n = (size_t)1 << logn;
	uint32_t *gm = tmp;
	uint32_t *t1 = gm + n;
	uint32_t *t2 = t1 + n;

	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	uint32_t R2 = PRIMES[0].R2;
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);

	/*
	 * t2 <- q + g*F (RNS+NTT)
	 */
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mp_set(g[u], p);
		t2[u] = mp_set(F[u], p);
	}
	mp_NTT(logn, t1, gm, p, p0i);
	mp_NTT(logn, t2, gm, p, p0i);
	uint32_t mq = mp_set(q, p);
	for (size_t u = 0; u < n; u ++) {
		uint32_t x = mp_montymul(t1[u], t2[u], p, p0i);
		t2[u] = mp_add(mq, mp_montymul(x, R2, p, p0i), p);
	}

	/*
	 * t2 <- (q + g*F)/f = G (RNS)
	 */
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mp_set(f[u], p);
	}
	mp_NTT(logn, t1, gm, p, p0i);
	uint32_t b = 0;
	for (size_t u = 0; u < n; u ++) {
		b |= t1[u] - 1;
		t2[u] = mp_div(t2[u], t1[u], p);
	}
	uint32_t *igm = gm;
	mp_mkigm(logn, igm, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t2, igm, p, p0i);

	/*
	 * Check that the G coefficients are in the proper range.
	 */
	int8_t *G = (int8_t *)tmp;
	for (size_t u = 0; u < n; u ++) {
		uint32_t x = t2[u];
		uint32_t y = tbmask((ulim << 1) - mp_add(x, ulim, p));
		b |= y;
		int32_t z = mp_norm(x & ~y, p);
		G[u] = (int8_t)z;
	}

	/*
	 * This failed if f was not invertible, i.e. one of its NTT
	 * coefficients was zero, of if any of the G coefficients was
	 * out of range.
	 */
	return (int)(1 - (b >> 31));
}

#if NTRUGEN_STATS
/*
 * All counters for statistics are gathered here. They are not thread-safe
 * and thus should be disabled in normal builds.
 */
uint32_t stats_hawk_ctt_attempt = 0;
uint32_t stats_hawk_ctt_reject = 0;
uint32_t stats_solve_attempt = 0;
uint32_t stats_solve_err_gcd = 0;
uint32_t stats_solve_err_reduce = 0;
uint32_t stats_solve_err_limit = 0;
uint32_t stats_solve_success = 0;
uint32_t stats_compute_w_attempt = 0;
uint32_t stats_compute_w_err_lim1 = 0;
uint32_t stats_compute_w_err_lim2 = 0;
uint32_t stats_compute_w_err_lim3 = 0;
uint32_t stats_compute_w_err_norm = 0;
uint32_t stats_compute_w_success = 0;

/* see ng_inner.h */
void
stats_init(void)
{
	stats_hawk_ctt_attempt = 0;
	stats_hawk_ctt_reject = 0;
	stats_solve_attempt = 0;
	stats_solve_err_gcd = 0;
	stats_solve_err_reduce = 0;
	stats_solve_err_limit = 0;
	stats_solve_success = 0;
	stats_compute_w_attempt = 0;
	stats_compute_w_err_lim1 = 0;
	stats_compute_w_err_lim2 = 0;
	stats_compute_w_err_lim3 = 0;
	stats_compute_w_err_norm = 0;
	stats_compute_w_success = 0;
}
#endif
