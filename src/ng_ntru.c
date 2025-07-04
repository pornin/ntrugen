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
 * Offset in tmp[] for saving the intermediate (f,g) values. It is
 * expressed in 32-bit words, and can use 'logn_top' to access the top-level
 * degree.
 */
#define FG_SAVE_OFFSET   ((size_t)5 << logn_top)

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
	if (!make_fg_deepest(prof, logn_top, f, g, tmp, FG_SAVE_OFFSET)) {
		return SOLVE_ERR_GCD;
	}

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
	 * (Only F is multiplied since G is ultimately discarded.)
	 */
	if (prof->q != 1) {
		if (zint_mul_small(Fp, len, prof->q) != 0) {
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
 * Input is F from one level deeper (half-degree), in plain
 * representation, at the start of tmp[]; output is F from this
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
	 *    input: hn * max_bl_small[depth + 1]
	 */

	unsigned logn = logn_top - depth;
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;

	/*
	 * slen   size for (f,g) at this level (also size of output F)
	 * llen   size for unreduced F at this level
	 * dlen   size for input F from deeper level
	 * Note: we always have llen >= dlen (constraint enforced in profiles)
	 */
	size_t slen = prof->max_bl_small[depth];
	size_t llen = prof->max_bl_large[depth];
	size_t dlen = prof->max_bl_small[depth + 1];

	/*
	 * Fd   F from deeper level (dlen*hn)
	 * ft   f from this level (slen*n)
	 */
	uint32_t *Fd = tmp;
	uint32_t *fgt = Fd + dlen * hn;

	/*
	 * Get (f,g) for this level (in RNS+NTT).
	 * Computation of the resultants left saved copies of (f,g) at
	 * the end of the tmp buffer, beyond some specific depth.
	 */
	if (depth < prof->min_save_fg[logn_top]) {
		make_fg_intermediate(prof, logn_top, f, g, depth, fgt);
	} else {
		uint32_t *sav_fg = tmp + FG_SAVE_OFFSET;
		for (unsigned d = prof->min_save_fg[logn_top];
			d <= depth; d ++)
		{
			sav_fg -= prof->max_bl_small[d] << (logn_top + 1 - d);
		}
		memmove(fgt, sav_fg, 2 * slen * n * sizeof *fgt);
	}

	/*
	 * Move buffers so that we have room for the unreduced (F,G) at
	 * this level.
	 *   Ft   F from this level (unreduced) (llen*n)
	 *   ft   f from this level (slen*n)
	 *   gt   g from this level (slen*n)
	 *   Fd   F from deeper level (dlen*hn)
	 */
	uint32_t *Ft = tmp;
	uint32_t *ft = Ft + llen * n;
	uint32_t *gt = ft + slen * n;
	Fd = gt + slen * n;
	uint32_t *t1 = Fd + dlen * hn;
	memmove(ft, fgt, 2 * n * slen * sizeof *ft);
	memmove(Fd, tmp, hn * dlen * sizeof *tmp);

	/*
	 * Convert Fd to RNS, with output temporarily stored in Ft. Fd
	 * has degree hn only; we store the values for each modulus p in
	 * the _last_ hn slots of the n-word line for that modulus.
	 */
	for (size_t u = 0; u < llen; u ++) {
		uint32_t p = PRIMES[u].p;
		uint32_t p0i = PRIMES[u].p0i;
		uint32_t R2 = PRIMES[u].R2;
		uint32_t Rx = mp_Rx31((unsigned)dlen, p, p0i, R2);
		uint32_t *xt = Ft + u * n + hn;
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
			}
		} else {
			for (size_t v = 0; v < hn; v ++) {
				xt[v] = zint_mod_small_signed(Fd + v, dlen, hn,
					p, p0i, R2, Rx);
			}
		}
#else // NTRUGEN_AVX2
		for (size_t v = 0; v < hn; v ++) {
			xt[v] = zint_mod_small_signed(Fd + v, dlen, hn,
				p, p0i, R2, Rx);
		}
#endif // NTRUGEN_AVX2
	}

	/*
	 * Fd is no longer needed.
	 */
	t1 = Fd;

	/*
	 * Compute F (unreduced) modulo sufficiently many small primes.
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
		 * Memory layout: we keep Ft, ft and gt; we append:
		 *   gm    NTT support (n)
		 *   igm   iNTT support (n)
		 *   gx    temporary g mod p (NTT) (n)
		 */
		uint32_t *gm = t1;
		uint32_t *igm = gm + n;
		uint32_t *gx = igm + n;
		mp_mkgmigm(logn, gm, igm, PRIMES[u].g, PRIMES[u].ig, p, p0i);
		if (u < slen) {
			memcpy(gx, gt + u * n, n * sizeof *gx);
			mp_iNTT(logn, ft + u * n, igm, p, p0i);
			mp_iNTT(logn, gt + u * n, igm, p, p0i);
		} else {
			uint32_t Rx = mp_Rx31((unsigned)slen, p, p0i, R2);
			for (size_t v = 0; v < n; v ++) {
				gx[v] = zint_mod_small_signed(gt + v, slen, n,
					p, p0i, R2, Rx);
			}
			mp_NTT(logn, gx, gm, p, p0i);
		}

		/*
		 * We have F from deeper level in Ft, in RNS. We apply the
		 * NTT modulo p.
		 */
		uint32_t *Fe = Ft + u * n;
		mp_NTT(logn - 1, Fe + hn, gm, p, p0i);

		/*
		 * Compute F (unreduced) modulo p.
		 */
#if NTRUGEN_AVX2
		if (hn >= 4) {
			__m256i yp = _mm256_set1_epi32(p);
			__m256i yp0i = _mm256_set1_epi32(p0i);
			__m256i yR2 = _mm256_set1_epi32(R2);
			for (size_t v = 0; v < hn; v += 4) {
				__m256i yga = _mm256_loadu_si256(
					(__m256i *)(gx + (v << 1)));
				__m256i ygb = _mm256_srli_epi64(yga, 32);
				__m128i xFe = _mm_loadu_si128(
					(__m128i *)(Fe + v + hn));
				__m256i yFp = _mm256_permute4x64_epi64(
					_mm256_castsi128_si256(xFe), 0x50);
				yFp = _mm256_shuffle_epi32(yFp, 0x30);
				yFp = mp_montymul_x4(yFp, yR2, yp, yp0i);
				__m256i yFe0 = mp_montymul_x4(
					ygb, yFp, yp, yp0i);
				__m256i yFe1 = mp_montymul_x4(
					yga, yFp, yp, yp0i);
				_mm256_storeu_si256((__m256i *)(Fe + (v << 1)),
					_mm256_or_si256(yFe0,
						_mm256_slli_epi64(yFe1, 32)));
			}
		} else {
			for (size_t v = 0; v < hn; v ++) {
				uint32_t ga = gx[(v << 1) + 0];
				uint32_t gb = gx[(v << 1) + 1];
				uint32_t mFp = mp_montymul(
					Fe[v + hn], R2, p, p0i);
				Fe[(v << 1) + 0] = mp_montymul(gb, mFp, p, p0i);
				Fe[(v << 1) + 1] = mp_montymul(ga, mFp, p, p0i);
			}
		}
#else // NTRUGEN_AVX2
		for (size_t v = 0; v < hn; v ++) {
			uint32_t ga = gx[(v << 1) + 0];
			uint32_t gb = gx[(v << 1) + 1];
			uint32_t mFp = mp_montymul(Fe[v + hn], R2, p, p0i);
			Fe[(v << 1) + 0] = mp_montymul(gb, mFp, p, p0i);
			Fe[(v << 1) + 1] = mp_montymul(ga, mFp, p, p0i);
		}
#endif // NTRUGEN_AVX2

		/*
		 * We want the new F in RNS only (no NTT).
		 */
		mp_iNTT(logn, Fe, igm, p, p0i);
	}

	/*
	 * We no longer need g.
	 */
	t1 = gt;

	/*
	 * Edge case: if slen == llen, then we have not rebuilt f
	 * into plain representation yet.
	 */
	if (slen == llen) {
		zint_rebuild_CRT(ft, slen, n, 1, 1, t1);
	}

	/*
	 * We now have the unreduced F in RNS. We rebuild its plain
	 * representation.
	 */
	zint_rebuild_CRT(Ft, llen, n, 1, 1, t1);

	/*
	 * We now reduce this F with Babai's round-off algorithm. The
	 * reduction conceptually goes as follows:
	 *   k <- round((F*adj(f) + G*adj(g))/(f*adj(f) + g*adj(g)))
	 *   (F, G) <- (F - k*f, G - k*g)
	 * We only have F; however, G is such that:
	 *   f*G - g*F = q
	 * hence:
	 *   G = (q + g*F)/f
	 * which we can move into the expression of k, and the expression
	 * simplifies into:
	 *   k = round(F/f + q*adj(g)/(f*(f*adj(f) + g*adj(g))))
	 * The second part only depends on f and g; moreover, it is
	 * heuristically negligible, i.e. we can compute an approximate
	 * value of k as:
	 *   k = round(F/f)
	 * which is a lot simpler. In practice the approximation is good
	 * enough for our purposes; a few coefficients are offset by +/-1
	 * but that does not prevent reduction from succeeding most of the
	 * time.
	 *
	 * We use fixed-point approximations of f and F to get a value k as
	 * a small polynomial with scaling; we then apply k on the
	 * full-width polynomial. Each iteration "shaves" a few bits off F.
	 *
	 * We apply the process sufficiently many times to reduce F
	 * to the size of f with a reasonable probability of success.
	 * Since we want full constant-time processing, the number of
	 * iterations and the accessed slots work on some assumptions on
	 * the sizes of values (sizes have been measured over many samples,
	 * and a margin of 5 times the standard deviation).
	 */

	/*
	 * If depth is at least 2, and we will use the NTT to subtract
	 * k*f from F, then we will need to convert f to NTT over
	 * slen+1 words, which requires an extra word to ft.
	 */
	int use_sub_ntt = (depth > 1 && logn >= MIN_LOGN_FGNTT);
	if (use_sub_ntt) {
		t1 += n;
	}

	/*
	 * New layout:
	 *   Ft    F from this level (unreduced) (llen*n)
	 *   ft    f from this level (slen*n) (+n if use_sub_ntt)
	 *   rt3   (n fxr = 2*n)
	 */
	fxr *rt3 = (fxr *)t1;

	/*
	 * We consider only the top rlen words of f.
	 */
	size_t rlen = prof->word_win[depth];
	if (rlen > slen) {
		rlen = slen;
	}
	size_t blen = slen - rlen;
	uint32_t *ftb = ft + blen * n;
	uint32_t scale_fg = 31 * (uint32_t)blen;
	uint32_t scale_FG = 31 * (uint32_t)llen;

	/*
	 * Convert f into fixed-point approximations, into rt3. It is scaled
	 * down by 2^(scale_fg + scale_x). scale_fg is public (it depends
	 * only on the recursion depth), but scale_x comes from a measurement
	 * on the actual coefficient values of f and is thus secret.
	 *
	 * The value scale_x is adjusted so that the largest coefficient is
	 * close to, but lower than, some limit t (in absolute value). The
	 * limit t is chosen so that f*adj(f) does not overflow, i.e. all
	 * coefficients must remain below 2^31.
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
	 * Since we want the coefficients of f*adj(f) not to exceed
	 * 2^31, we need n^2*t^2 <= 2^31, i.e. n*t <= 2^15.5. We can adjust t
	 * accordingly (called scale_t in the code below). We also need to
	 * take care that t must not exceed scale_x. Approximation of f
	 * is extracted with scale scale_fg + scale_x - scale_t, and
	 * later fixed by dividing it by 2^scale_t.
	 */
	uint32_t scale_x = poly_max_bitlength(logn, ftb, rlen);
	uint32_t scale_t = 15 - logn;
	scale_t ^= (scale_t ^ scale_x) & tbmask(scale_x - scale_t);
	uint32_t scdiff = scale_x - scale_t;

	poly_big_to_fixed(logn, rt3, ftb, rlen, scdiff);

	/*
	 * Compute adj(f)/(f*adj(f)) into rt3 (FFT, with scaling).
	 */
	vect_FFT(logn, rt3);
	vect_inv_mul2e_fft(logn, rt3, scale_t);

	/*
	 * New layout:
	 *   Ft    F from this level (unreduced) (llen*n)
	 *   ft    f from this level (slen*n) (+n if use_sub_ntt)
	 *   rt3   (n fxr = 2*n)
	 *   rt1   (n fxr = 2*n)     |   k    (n)
	 *                           |   t2   (3*n)
	 */
	fxr *rt1 = rt3 + n;
	int32_t *k = (int32_t *)rt1;
	uint32_t *t2 = (uint32_t *)(k + n);

	/*
	 * If we are going to use poly_sub_scaled_ntt(), then we convert
	 * f to the NTT representation. Since poly_sub_scaled_ntt()
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
	}

	/*
	 * Reduce F repeatedly.
	 */
	size_t FGlen = llen;
	for (;;) {
		/*
		 * Convert the current F into fixed-point. We want
		 * to apply scaling scale_FG + scale_x.
		 */
		uint32_t tlen, toff;
		DIVREM31(tlen, toff, scale_FG);
		poly_big_to_fixed(logn, rt1,
			Ft + tlen * n, FGlen - tlen, scale_x + toff);

		/*
		 * rt1 <- (F*adj(f)) / (f*adj(f))
		 */
		vect_FFT(logn, rt1);
		vect_mul_fft(logn, rt1, rt3);
		vect_iFFT(logn, rt1);

		/*
		 * k <- round(rt1)
		 */
		for (size_t u = 0; u < n; u ++) {
			k[u] = fxr_round(rt1[u]);
		}

		/*
		 * f is scaled by scale_fg + scale_x
		 * F is scaled by scale_FG + scale_x
		 * Thus, k is scaled by scale_FG - scale_fg, which is public.
		 */
		uint32_t scale_k = scale_FG - scale_fg;
		if (depth == 1) {
			poly_sub_kf_scaled_depth1(logn_top, Ft, FGlen,
				(uint32_t *)k, scale_k, f, t2);
		} else if (use_sub_ntt) {
			poly_sub_scaled_ntt(logn, Ft, FGlen, ft, slen,
				k, scale_k, t2);
		} else {
			poly_sub_scaled(logn, Ft, FGlen, ft, slen, k, scale_k);
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
			/*
			 * We decrement FGlen; when we do so, we check that
			 * it does not damage any of the value, i.e. that the
			 * removed words are redundant with the the remaining
			 * words. If a non-redundant word is removed then
			 * the reduction failed. A contrario, as long as such
			 * removal works, then we know that the reduction
			 * keeps working and the computation so far is
			 * correct.
			 */
			FGlen --;
			uint32_t *xp = &Ft[(FGlen - 1) << logn];
			for (size_t u = 0; u < n; u ++) {
				uint32_t sw = -(xp[u] >> 30) >> 1;
				if (xp[u + n] != sw) {
					return SOLVE_ERR_REDUCE;
				}
			}
		}
	}

	/*
	 * Output F is already in the right place.
	 */
	return SOLVE_OK;
}

/*
 * Solving the NTRU equation, top recursion level. This is a specialized
 * variant for solve_NTRU_intermediate() with depth == 0, for lower RAM
 * usage and faster operation; it also ensures better precision for the
 * reduction. This function returns both F and G, in that order, at the
 * start of tmp[].
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
	 * On input, Fd from upper level (hn words) is at the start of
	 * tmp[].
	 */
	uint32_t *t1 = tmp;
	uint32_t *t2 = t1 + n;
	uint32_t *t3 = t2 + n;
	uint32_t *t4 = t3 + n;
	uint32_t *t5 = t4 + n;

	/*
	 * Convert Fd to RNS+NTT, into t3.
	 */
	uint32_t *gm = t4;
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	poly_mp_set(logn - 1, t1, p);
	mp_NTT(logn - 1, t1, gm, p, p0i);
	memcpy(t3, t1, hn * sizeof *t1);

	/*
	 * Compute F (unreduced, RNS+NTT) into t1.
	 */
	poly_mp_set_small(logn, t2, g, p);
	mp_NTT(logn, t2, gm, p, p0i);
#if NTRUGEN_AVX2
	if (hn >= 4) {
		__m256i yp = _mm256_set1_epi32(p);
		__m256i yp0i = _mm256_set1_epi32(p0i);
		__m256i yR2 = _mm256_set1_epi32(R2);
		for (size_t u = 0; u < hn; u += 4) {
			__m256i yga = _mm256_loadu_si256(
				(__m256i *)(t2 + (u << 1)));
			__m256i ygb = _mm256_srli_epi64(yga, 32);
			__m128i xFd = _mm_loadu_si128((__m128i *)(t3 + u));
			__m256i yFd = _mm256_permute4x64_epi64(
				_mm256_castsi128_si256(xFd), 0x50);
			yFd = _mm256_shuffle_epi32(yFd, 0x30);
			yFd = mp_montymul_x4(yFd, yR2, yp, yp0i);
			__m256i yFe0 = mp_montymul_x4(ygb, yFd, yp, yp0i);
			__m256i yFe1 = mp_montymul_x4(yga, yFd, yp, yp0i);
			_mm256_storeu_si256((__m256i *)(t1 + (u << 1)),
				_mm256_or_si256(yFe0,
					_mm256_slli_epi64(yFe1, 32)));
		}
	} else {
		for (size_t u = 0; u < hn; u ++) {
			uint32_t ga = t2[(u << 1) + 0];
			uint32_t gb = t2[(u << 1) + 1];
			uint32_t mF = mp_montymul(t3[u], R2, p, p0i);
			t1[(u << 1) + 0] = mp_montymul(gb, mF, p, p0i);
			t1[(u << 1) + 1] = mp_montymul(ga, mF, p, p0i);
		}
	}
#else // NTRUGEN_AVX2
	for (size_t u = 0; u < hn; u ++) {
		uint32_t ga = t2[(u << 1) + 0];
		uint32_t gb = t2[(u << 1) + 1];
		uint32_t mF = mp_montymul(t3[u], R2, p, p0i);
		t1[(u << 1) + 0] = mp_montymul(gb, mF, p, p0i);
		t1[(u << 1) + 1] = mp_montymul(ga, mF, p, p0i);
	}
#endif

	/*
	 * Layout:
	 *   t1   F (unreduced, RNS+NTT)
	 *   t2   g (RNS+NTT)
	 *   t3   free
	 *   t4   gm (NTT support)
	 *   t5   free
	 */

	/*
	 * Convert f to RNS+NTT. Since we are going to divide by f modulo p,
	 * we also need to check that f is invertible modulo p (which should
	 * almost always be the case in practice).
	 *   t3 <- f (RNS+NTT)
	 */
	poly_mp_set_small(logn, t3, f, p);
	mp_NTT(logn, t3, gm, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		if (t3[u] == 0) {
			return SOLVE_ERR_REDUCE;
		}
	}

	/*
	 * Layout:
	 *   t1   F (unreduced, RNS+NTT)
	 *   t2   g (RNS+NTT)
	 *   t3   f (RNS+NTT)
	 *   t4   free
	 *   t5   free
	 */

	/*
	 * We want to perform the reduction. Since this is the last one,
	 * we want to be precise, i.e. to use the full expression for k:
	 *
	 *   k = round((F*adj(f) + G*adj(g))/(f*adj(f) + g*adj(g)))
	 *
	 * We do not have G but we know that G = (q + g*F)/f, which we
	 * can compute modulo p (the division by f is exact over the
	 * integers, hence computing it modulo p yields the correct result,
	 * as long as the coefficients of G are in [-p/2,+p/2], which is
	 * heuristically the case). We accumulate the numerator and
	 * denominator into t2 and t3, respectively.
	 */
	uint32_t q = prof->q;
#if NTRUGEN_AVX2
	if (hn >= 8) {
		__m256i yp = _mm256_set1_epi32(p);
		__m256i yp0i = _mm256_set1_epi32(p0i);
		__m256i yR2 = _mm256_set1_epi32(R2);
		__m256i yq = _mm256_set1_epi32(q);
		for (size_t u = 0; u < hn; u += 8) {
			/* Load values by groups of 8. */
			__m256i yf0 = _mm256_loadu_si256(
				(__m256i *)(t3 + u));
			__m256i yf1 = _mm256_loadu_si256(
				(__m256i *)(t3 + (n - 8 - u)));
			__m256i yg0 = _mm256_loadu_si256(
				(__m256i *)(t2 + u));
			__m256i yg1 = _mm256_loadu_si256(
				(__m256i *)(t2 + (n - 8 - u)));
			__m256i yF0 = _mm256_loadu_si256(
				(__m256i *)(t1 + u));
			__m256i yF1 = _mm256_loadu_si256(
				(__m256i *)(t1 + (n - 8 - u)));
			__m256i ymg0 = mp_montymul_x8(yg0, yR2, yp, yp0i);
			__m256i ymg1 = mp_montymul_x8(yg1, yR2, yp, yp0i);

			/* Compute G. */
			__m256i yG0 = mp_div_x8(
				mp_add_x8(mp_montymul_x8(ymg0, yF0, yp, yp0i),
					yq, yp),
				yf0, yp);
			__m256i yG1 = mp_div_x8(
				mp_add_x8(mp_montymul_x8(ymg1, yF1, yp, yp0i),
					yq, yp),
				yf1, yp);

			/* We need adj(f) and adj(g), which entails reversing
			   the order of the values. */
			__m256i yfa1 = _mm256_permute4x64_epi64(yf0, 0x4E);
			yfa1 = _mm256_shuffle_epi32(yfa1, 0x1B);
			__m256i yfa0 = _mm256_permute4x64_epi64(yf1, 0x4E);
			yfa0 = _mm256_shuffle_epi32(yfa0, 0x1B);
			__m256i ymga1 = _mm256_permute4x64_epi64(ymg0, 0x4E);
			ymga1 = _mm256_shuffle_epi32(ymga1, 0x1B);
			__m256i ymga0 = _mm256_permute4x64_epi64(ymg1, 0x4E);
			ymga0 = _mm256_shuffle_epi32(ymga0, 0x1B);

			__m256i ymfa0 = mp_montymul_x8(yfa0, yR2, yp, yp0i);
			__m256i ymfa1 = mp_montymul_x8(yfa1, yR2, yp, yp0i);

			/* kn <- F*adj(f) + G*adj(g) */
			__m256i ykn0 = mp_add_x8(
				mp_montymul_x8(ymfa0, yF0, yp, yp0i),
				mp_montymul_x8(ymga0, yG0, yp, yp0i), yp);
			__m256i ykn1 = mp_add_x8(
				mp_montymul_x8(ymfa1, yF1, yp, yp0i),
				mp_montymul_x8(ymga1, yG1, yp, yp0i), yp);

			/* kd <- f*adj(f) + g*adj(g)
			   kd is auto-adjoint, hence we can compute the
			   right half as a swap of the left half. */
			__m256i ykd0 = mp_add_x8(
				mp_montymul_x8(ymfa0, yf0, yp, yp0i),
				mp_montymul_x8(ymga0, yg0, yp, yp0i), yp);
			__m256i ykd1 = _mm256_permute4x64_epi64(ykd0, 0x4E);
			ykd1 = _mm256_shuffle_epi32(ykd1, 0x1B);

			_mm256_storeu_si256((__m256i *)(t2 + u), ykn0);
			_mm256_storeu_si256((__m256i *)(t2 + n - 8 - u), ykn1);
			_mm256_storeu_si256((__m256i *)(t3 + u), ykd0);
			_mm256_storeu_si256((__m256i *)(t3 + n - 8 - u), ykd1);
		}
	} else {
		for (size_t u = 0; u < hn; u ++) {
			uint32_t tf0 = t3[u];
			uint32_t tf1 = t3[n - 1 - u];
			uint32_t tg0 = t2[u];
			uint32_t tg1 = t2[n - 1 - u];
			uint32_t tF0 = t1[u];
			uint32_t tF1 = t1[n - 1 - u];
			uint32_t mf0 = mp_montymul(tf0, R2, p, p0i);
			uint32_t mf1 = mp_montymul(tf1, R2, p, p0i);
			uint32_t mg0 = mp_montymul(tg0, R2, p, p0i);
			uint32_t mg1 = mp_montymul(tg1, R2, p, p0i);
			uint32_t tG0 = mp_div(
				mp_add(q, mp_montymul(mg0, tF0, p, p0i), p),
				tf0, p);
			uint32_t tG1 = mp_div(
				mp_add(q, mp_montymul(mg1, tF1, p, p0i), p),
				tf1, p);
			uint32_t kn0 = mp_add(
				mp_montymul(mf1, tF0, p, p0i),
				mp_montymul(mg1, tG0, p, p0i), p);
			uint32_t kn1 = mp_add(
				mp_montymul(mf0, tF1, p, p0i),
				mp_montymul(mg0, tG1, p, p0i), p);
			uint32_t kd = mp_add(
				mp_montymul(mf0, tf1, p, p0i),
				mp_montymul(mg0, tg1, p, p0i), p);
			t2[u] = kn0;
			t2[n - 1 - u] = kn1;
			t3[u] = kd;
			t3[n - 1 - u] = kd;
		}
	}
#else // NTRUGEN_AVX2
	for (size_t u = 0; u < hn; u ++) {
		uint32_t tf0 = t3[u];
		uint32_t tf1 = t3[n - 1 - u];
		uint32_t tg0 = t2[u];
		uint32_t tg1 = t2[n - 1 - u];
		uint32_t tF0 = t1[u];
		uint32_t tF1 = t1[n - 1 - u];
		uint32_t mf0 = mp_montymul(tf0, R2, p, p0i);
		uint32_t mf1 = mp_montymul(tf1, R2, p, p0i);
		uint32_t mg0 = mp_montymul(tg0, R2, p, p0i);
		uint32_t mg1 = mp_montymul(tg1, R2, p, p0i);
		uint32_t tG0 = mp_div(
			mp_add(q, mp_montymul(mg0, tF0, p, p0i), p), tf0, p);
		uint32_t tG1 = mp_div(
			mp_add(q, mp_montymul(mg1, tF1, p, p0i), p), tf1, p);
		uint32_t kn0 = mp_add(
			mp_montymul(mf1, tF0, p, p0i),
			mp_montymul(mg1, tG0, p, p0i), p);
		uint32_t kn1 = mp_add(
			mp_montymul(mf0, tF1, p, p0i),
			mp_montymul(mg0, tG1, p, p0i), p);
		uint32_t kd = mp_add(
			mp_montymul(mf0, tf1, p, p0i),
			mp_montymul(mg0, tg1, p, p0i), p);
		t2[u] = kn0;
		t2[n - 1 - u] = kn1;
		t3[u] = kd;
		t3[n - 1 - u] = kd;
	}
#endif

	/*
	 * Current layout:
	 *   t1   unreduced F (RNS+NTT)
	 *   t2   F*adj(f) + G*adj(g) (RNS+NTT)
	 *   t3   f*adj(f) + g*adj(g) (RNS+NTT)
	 *   t4   free
	 *   t5   free
	 */

	/*
	 * Convert back numerator and denominator to plain integers.
	 */
	mp_mkigm(logn, t4, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t2, t4, p, p0i);
	mp_iNTT(logn, t3, t4, p, p0i);
#if NTRUGEN_AVX2
	if (n >= 8) {
		__m256i yp = _mm256_set1_epi32(p);
		__m256i yhp = _mm256_set1_epi32((p + 1) >> 1);
		for (size_t u = 0; u < n; u += 8) {
			__m256i yn = _mm256_loadu_si256((__m256i *)(t2 + u));
			__m256i yd = _mm256_loadu_si256((__m256i *)(t3 + u));
			yn = mp_norm_x8(yn, yp, yhp);
			yd = mp_norm_x8(yd, yp, yhp);
			_mm256_storeu_si256((__m256i *)(t2 + u), yn);
			_mm256_storeu_si256((__m256i *)(t3 + u), yd);
		}
	} else {
		for (size_t u = 0; u < n; u ++) {
			t2[u] = (uint32_t)mp_norm(t2[u], p);
			t3[u] = (uint32_t)mp_norm(t3[u], p);
		}
	}
#else // NTRUGEN_AVX2
	for (size_t u = 0; u < n; u ++) {
		/*
		 * NOTE: no truncature to 31 bits.
		 */
		t2[u] = (uint32_t)mp_norm(t2[u], p);
		t3[u] = (uint32_t)mp_norm(t3[u], p);
	}
#endif

#define DOWNSCALE   10
	/*
	 * We need to divide t2 by t3, and round the result. We convert
	 * them to FFT representation, downscaled by 2^10 (to avoid overflows).
	 * We first convert f*adj(f) + g*adj(g), which is auto-adjoint;
	 * thus, its FFT representation only has half-size.
	 */
	fxr *rt4 = (fxr *)t4;
	for (size_t u = 0; u < n; u ++) {
		uint64_t x = (uint64_t)*(int32_t *)&t3[u] << (32 - DOWNSCALE);
		rt4[u] = fxr_of_scaled32(x);
	}
	vect_FFT(logn, rt4);
	memcpy(rt4 + hn, rt4, hn * sizeof *rt4);
	fxr *rt5 = rt4 + hn;
	fxr *rt3 = (fxr *)t3;
	for (size_t u = 0; u < n; u ++) {
		uint64_t x = (uint64_t)*(int32_t *)&t2[u] << (32 - DOWNSCALE);
		rt3[u] = fxr_of_scaled32(x);
	}
	vect_FFT(logn, rt3);
#undef DOWNSCALE

	/*
	 * Current layout:
	 *   t1   unreduced F (RNS+NTT)
	 *   t2   free
	 *   t3   F*adj(f) + G*adj(g) (FFT) (first half)   <- alias: rt3
	 *   t4   F*adj(f) + G*adj(g) (FFT) (second half)
	 *   t5   f*adj(f) + g*adj(g) (FFT) (half-size)    <- alias: rt5
	 */

	/*
	 * Divide F*adj(f) + G*adj(g) by f*adj(f) + g*adj(g) and round
	 * the result into t2, with conversion to RNS.
	 */
	vect_div_autoadj_fft(logn, rt3, rt5);
	vect_iFFT(logn, rt3);
	for (size_t u = 0; u < n; u ++) {
		t2[u] = mp_set(fxr_round(rt3[u]), p);
	}

	/*
	 * Current layout:
	 *   t1   unreduced F (RNS+NTT)
	 *   t2   k (RNS)
	 *   t3   free
	 *   t4   free
	 *   t5   free
	 */

	/*
	 * Get back f and g, and convert all polynomials to RNS+NTT.
	 */
	gm = t5;
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	poly_mp_set_small(logn, t3, f, p);
	poly_mp_set_small(logn, t4, g, p);
	mp_NTT(logn, t2, gm, p, p0i);
	mp_NTT(logn, t3, gm, p, p0i);
	mp_NTT(logn, t4, gm, p, p0i);

	/*
	 * Current layout:
	 *   t1   unreduced F (RNS+NTT)
	 *   t2   k (RNS+NTT)
	 *   t3   f (RNS+NTT)
	 *   t4   g (RNS+NTT)
	 *   t5   free
	 */

	/*
	 * Reduce F by subtracting k*f, and recompute the corresponding G
	 * with:
	 *   G = (q + g*F)/f
	 * (We computed the unreduced G, but did not keep it, in order to
	 * save RAM; if we had kept it in an extra buffer, then we could
	 * reduce it by subtracting k*g, which would be faster than the
	 * division.)
	 */
#if NTRUGEN_AVX2
	if (n >= 8) {
		__m256i yp = _mm256_set1_epi32(p);
		__m256i yp0i = _mm256_set1_epi32(p0i);
		__m256i yR2 = _mm256_set1_epi32(R2);
		__m256i yq = _mm256_set1_epi32(q);
		for (size_t u = 0; u < n; u += 8) {
			__m256i yF = _mm256_loadu_si256((__m256i *)(t1 + u));
			__m256i yk = _mm256_loadu_si256((__m256i *)(t2 + u));
			__m256i yf = _mm256_loadu_si256((__m256i *)(t3 + u));
			__m256i yg = _mm256_loadu_si256((__m256i *)(t4 + u));
			__m256i ymf = mp_montymul_x8(yf, yR2, yp, yp0i);
			__m256i ymg = mp_montymul_x8(yg, yR2, yp, yp0i);
			yF = mp_sub_x8(yF,
				mp_montymul_x8(ymf, yk, yp, yp0i), yp);
			__m256i yG = mp_div_x8(
				mp_add_x8(yq,
					mp_montymul_x8(ymg, yF, yp, yp0i), yp),
				yf, yp);
			_mm256_storeu_si256((__m256i *)(t1 + u), yF);
			_mm256_storeu_si256((__m256i *)(t2 + u), yG);
		}
	} else {
		for (size_t u = 0; u < n; u ++) {
			uint32_t tF = t1[u];
			uint32_t tk = t2[u];
			uint32_t tf = t3[u];
			uint32_t tg = t4[u];
			uint32_t mf = mp_montymul(tf, R2, p, p0i);
			uint32_t mg = mp_montymul(tg, R2, p, p0i);
			tF = mp_sub(tF, mp_montymul(mf, tk, p, p0i), p);
			uint32_t tG = mp_div(
				mp_add(q, mp_montymul(mg, tF, p, p0i), p),
				tf, p);
			t1[u] = tF;
			t2[u] = tG;
		}
	}
#else // NTRUGEN_AVX2
	for (size_t u = 0; u < n; u ++) {
		uint32_t tF = t1[u];
		uint32_t tk = t2[u];
		uint32_t tf = t3[u];
		uint32_t tg = t4[u];
		uint32_t mf = mp_montymul(tf, R2, p, p0i);
		uint32_t mg = mp_montymul(tg, R2, p, p0i);
		tF = mp_sub(tF, mp_montymul(mf, tk, p, p0i), p);
		uint32_t tG = mp_div(
			mp_add(q, mp_montymul(mg, tF, p, p0i), p), tf, p);
		t1[u] = tF;
		t2[u] = tG;
	}
#endif

	/*
	 * Convert back F and G into normal representation.
	 */
	mp_mkigm(logn, t3, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t1, t3, p, p0i);
	mp_iNTT(logn, t2, t3, p, p0i);
	poly_mp_norm(logn, t1, p);
	poly_mp_norm(logn, t2, p);

	/*
	 * We're done! By construction, f*G - g*F = q modulo p; if both
	 * F and G are in the correct range ([-127,+127]) then this
	 * equation will also hold over plain integers:
	 *   N_inf(f*G - g*F) <= (127^2)*n*2 < 2^25 < p/2
	 * Verifying that F and G are in range is done by the caller.
	 */
	return SOLVE_OK;
}

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
