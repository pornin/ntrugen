#include "ng_inner.h"

const ntru_profile SOLVE_BAT_128_256 = {
	128,
	8, 8,
	{ 1, 1, 1, 2, 3, 4, 8, 14, 27, 0, 0 },
	{ 1, 1, 2, 3, 6, 11, 20, 39, 0, 0 },
	{ 1, 1, 1, 2, 3, 3, 3, 4, 0, 0 },
	13,
	{ 0, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};

const ntru_profile SOLVE_BAT_257_512 = {
	257,
	9, 9,
	{ 1, 1, 1, 2, 3, 5, 8, 15, 30, 59, 0 },
	{ 1, 1, 2, 4, 6, 12, 23, 44, 87, 0 },
	{ 1, 1, 1, 2, 3, 3, 3, 4, 6, 0 },
	13,
	{ 0, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};

const ntru_profile SOLVE_BAT_769_1024 = {
	769,
	10, 10,
	{ 1, 1, 1, 2, 3, 5, 10, 18, 35, 69, 137 },
	{ 1, 1, 2, 4, 7, 14, 27, 52, 102, 204 },
	{ 1, 1, 1, 2, 3, 3, 3, 4, 5, 7 },
	11,
	{ 0, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};

/* BAT, q = 128, n = 256 -> kmax = 2 */
const uint16_t gauss_BAT_128_256[] = {
	2,
	  153, 10823, 54712, 65382
};

/* BAT, q = 257, n = 512 -> kmax = 2 */
const uint16_t gauss_BAT_257_512[] = {
	2,
	  157, 10865, 54670, 65378
};

/* BAT, q = 769, n = 1024 -> kmax = 3 */
const uint16_t gauss_BAT_769_1024[] = {
	3,
	    1,   400, 12951, 52584, 65135, 65534
};

/*
 * Compute the vector w:
 *   w = round(q'*(gamma^2*F*adj(f) + G*adj(g))/(gamma^2*f*adj(f) + g*adj(g)))
 * with:
 *   q' = 64513
 *   gamma = 1 (for logn = 8 or 9) or sqrt(5) (for logn = 10)
 * Returned value: 1 on success, or 0 on error.
 *
 * Value of w is written at the start of tmp[] (in plain 32-bit format).
 *
 * Errors include an invalid parameter (logn must be 8, 9 or 10), or a
 * computation overflow. In the latter case, the function may return early;
 * thus, for constant-time discipline, any error should induce rejection
 * of the (f,g,F,G) quadruplet.
 *
 * RAM USAGE: 4.5*n words
 */
static int
compute_w(unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	uint32_t *restrict tmp)
{
	uint32_t qp = 64513;

	int32_t gamma2;
	uint64_t max_dnorm;
	switch (logn) {
	case 8:
		gamma2 = 1;
		max_dnorm = 756381852074;
		break;
	case 9:
		gamma2 = 1;
		max_dnorm = 1524693345421;
		break;
	case 10:
		gamma2 = 5;
		max_dnorm = 11296783533487;
		break;
	default:
		return 0;
	}

	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
#if NTRUGEN_STATS
	stats_compute_w_attempt ++;
#endif

	/*
	 * We create four buffers of size n words and one extra buffer of
	 * size n/2 words.
	 */
	uint32_t *t1 = tmp;
	uint32_t *t2 = t1 + n;
	uint32_t *t3 = t2 + n;
	uint32_t *t4 = t3 + n;
	uint32_t *t5 = t4 + n;   /* half-size */

	/*
	 * We want:
	 *   t1 <- (gamma^2)*F*adj(f) + G*adj(g)  (RNS+NTT)
	 *   t2 <- (gamma^2)*f*adj(f) + g*adj(g)  (RNS+NTT)
	 * We work modulo a big 31-bit prime, which is large enough to
	 * avoid wrap-arounds, and thus allows us to get the plain integer
	 * values (after normalization).
	 */
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	uint32_t R2 = PRIMES[0].R2;
	uint32_t *gm = t1;
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);

	/*
	 * gmv <- R*gamma^2  (Montgomery representation of gamma^2)
	 */
	uint32_t gmv = mp_montymul(R2, mp_set(gamma2, p), p, p0i);

	/*
	 * t2 <- f  (RNS+NTT)
	 * t3 <- g  (RNS+NTT)
	 * t4 <- F  (RNS+NTT)
	 */
	poly_mp_set_small(logn, t2, f, p);
	poly_mp_set_small(logn, t3, g, p);
	poly_mp_set_small(logn, t4, F, p);
	mp_NTT(logn, t2, gm, p, p0i);
	mp_NTT(logn, t3, gm, p, p0i);
	mp_NTT(logn, t4, gm, p, p0i);

	/*
	 * t5 <- (gamma^2)*f*adj(f) + g*adj(g)  (RNS+NTT)  (auto-adjoint)
	 * t4 <- (gamma^2)*F*adj(f)             (RNS+NTT)
	 */
	for (size_t u = 0; u < hn; u ++) {
		uint32_t xf = t2[u];
		uint32_t xfa = t2[n - 1 - u];
		uint32_t xg = t3[u];
		uint32_t xga = t3[n - 1 - u];
		uint32_t xF = t4[u];
		uint32_t xFa = t4[n - 1 - u];
		uint32_t gmvf = mp_montymul(gmv, xf, p, p0i);
		uint32_t gmvfa = mp_montymul(gmv, xfa, p, p0i);
		t5[u] = mp_montymul(R2, mp_add(
			mp_montymul(xf, gmvfa, p, p0i),
			mp_montymul(xg, xga, p, p0i), p), p, p0i);
		t4[u] = mp_montymul(R2,
			mp_montymul(xF, gmvfa, p, p0i), p, p0i);
		t4[n - 1 - u] = mp_montymul(R2,
			mp_montymul(xFa, gmvf, p, p0i), p, p0i);
	}

	/*
	 * t1 <- (gamma^2)*F*adj(f) + G*adj(g)  (RNS+NTT)
	 */
	poly_mp_set_small(logn, t2, G, p);
	mp_NTT(logn, t2, gm, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		uint32_t xG = t2[u];
		uint32_t xga = t3[n - 1 - u];
		uint32_t xz = t4[u];
		uint32_t mga = mp_montymul(xga, R2, p, p0i);
		t1[u] = mp_add(mp_montymul(mga, xG, p, p0i), xz, p);
	}

	/*
	 * t2 <- (gamma^2)*f*adj(f) + g*adj(j)  (RNS+NTT) (full-size)
	 */
	for (size_t u = 0; u < hn; u ++) {
		t2[u] = t2[n - 1 - u] = t5[u];
	}

	/*
	 * Convert t1 and t2 to plain (32-bit).
	 */
	uint32_t *igm = t3;
	mp_mkigm(logn, igm, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t1, igm, p, p0i);
	mp_iNTT(logn, t2, igm, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = (uint32_t)mp_norm(t1[u], p);
		t2[u] = (uint32_t)mp_norm(t2[u], p);
	}

	/*
	 * For the division, we go into the FFT domain. We check that the
	 * FFT won't overflow. The dividend is scaled down by 10 bits to
	 * compensate for the multiplication by q'.
	 *
	 * Buffer reorganization:
	 *    t1    (gamma^2)*F*adj(f) + G*adj(g) (plain, 32-bit) (n)
	 *    t2    (gamma^2)*f*adj(f) + g*adj(g) (plain, 32-bit) (n)
	 *    rt1   receives the dividend (FFT) (n fxr = 2*n)
	 *
	 * rt2 is an alias on t2 (as hn fxr values) and recieves the
	 * divisor (in FFT). The divisor is auto-adjoint, so it uses
	 * only half of the space in FFT representation; we compute in
	 * rt1 then move it to rt2.
	 */
	fxr *rt1 = (fxr *)(t2 + n);
	for (size_t u = 0; u < n; u ++) {
		rt1[u] = fxr_of(*(int32_t *)&t2[u]);
	}
	vect_FFT(logn, rt1);
	fxr *rt2 = (fxr *)t2;
	memmove(rt2, rt1, hn * sizeof *rt1);

	/*
	 * For the dividend, we multiply by q' but also scale down
	 * by 2^10; we check that the operation won't overflow.
	 */
	int32_t lim1 = (int32_t)(((uint64_t)1 << (41 - logn)) / qp);
	for (size_t u = 0; u < n; u ++) {
		int32_t x = *(int32_t *)&t1[u];
		if (x <= -lim1 || x >= +lim1) {
#if NTRUGEN_STATS
			stats_compute_w_err_lim1 ++;
#endif
			return 0;
		}
		rt1[u] = fxr_of_scaled32(((uint64_t)x * qp) << 22);
	}
	vect_FFT(logn, rt1);

	/*
	 * Divisor is auto-adjoint. We inline the division loop here because
	 * we also want to check on overflows.
	 */
	for (size_t u = 0; u < hn; u ++) {
		fxr z1r = rt1[u];
		fxr z1i = rt1[u + hn];
		fxr z2 = rt2[u];
		if (!fxr_lt(fxr_div2e(fxr_abs(z1r), 30 - logn), z2)
			|| !fxr_lt(fxr_div2e(fxr_abs(z1i), 30 - logn), z2))
		{
#if NTRUGEN_STATS
			stats_compute_w_err_lim2 ++;
#endif
			return 0;
		}
		rt1[u] = fxr_div(z1r, z2);
		rt1[u + hn] = fxr_div(z1i, z2);
	}
	vect_iFFT(logn, rt1);

	/*
	 * The unrounded w is in rt1 (scaled down by 2^10); we just have to
	 * round the coefficients and check that they are all in the
	 * allowed [-2^16..+2^16] range.
	 */
	fxr lim2 = fxr_of(1 << 6);
	for (size_t u = 0; u < n; u ++) {
		if (fxr_lt(lim2, fxr_abs(rt1[u]))) {
#if NTRUGEN_STATS
			stats_compute_w_err_lim3 ++;
#endif
			return 0;
		}
		t1[u] = (uint32_t)fxr_round(fxr_mul2e(rt1[u], 10));
	}

	/*
	 * Check that the norm of (gamma*Fd, Gd) is low enough, with:
	 *   Fd = q'*F - f*w
	 *   Gd = q'*G - G*w
	 *
	 * Buffer layout:
	 *   t1    w (plain, 32-bit) (n)
	 *   t2    free (n)
	 *   t3    free (n)
	 *   t4    free (n)
	 * gm and igm are set to both point to t4 (for space-saving reasons,
	 * we use the same buffer for both values).
	 */
	gm = igm = t4;

	/*
	 * Convert w to NTT + Montgomery.
	 */
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mp_montymul(R2, mp_set(*(int32_t *)&t1[u], p), p, p0i);
	}
	mp_NTT(logn, t1, gm, p, p0i);

	/*
	 * t2 <- f     (NTT)
	 * t3 <- q'*F  (NTT)
	 */
	for (size_t u = 0; u < n; u ++) {
		t2[u] = mp_set(f[u], p);
		t3[u] = mp_set((int32_t)qp * F[u], p);
	}
	mp_NTT(logn, t2, gm, p, p0i);
	mp_NTT(logn, t3, gm, p, p0i);

	/*
	 * t2 <- Fd, and compute its squared norm.
	 */
	for (size_t u = 0; u < n; u ++) {
		t2[u] = mp_sub(t3[u], mp_montymul(t2[u], t1[u], p, p0i), p);
	}
	mp_mkigm(logn, igm, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t2, igm, p, p0i);
	uint64_t Fdnorm = 0;
	for (size_t u = 0; u < n; u ++) {
		int32_t x = mp_norm(t2[u], p);
		Fdnorm += (uint64_t)((int64_t)x * (int64_t)x);
	}

	/*
	 * t2 <- g     (NTT)
	 * t3 <- q'*G  (NTT)
	 */
	for (size_t u = 0; u < n; u ++) {
		t2[u] = mp_set(g[u], p);
		t3[u] = mp_set((int32_t)qp * G[u], p);
	}
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	mp_NTT(logn, t2, gm, p, p0i);
	mp_NTT(logn, t3, gm, p, p0i);

	/*
	 * t2 <- Gd, and compute its squared norm.
	 */
	for (size_t u = 0; u < n; u ++) {
		t2[u] = mp_sub(t3[u], mp_montymul(t2[u], t1[u], p, p0i), p);
	}
	mp_mkigm(logn, igm, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t2, igm, p, p0i);
	uint64_t Gdnorm = 0;
	for (size_t u = 0; u < n; u ++) {
		int32_t x = mp_norm(t2[u], p);
		Gdnorm += (uint64_t)((int64_t)x * (int64_t)x);
	}

	/*
	 * Check that the total squared norm of (Fd, gamma*Gd) is low enough.
	 */
	if ((uint64_t)gamma2 * Fdnorm + Gdnorm > max_dnorm) {
#if NTRUGEN_STATS
		stats_compute_w_err_norm ++;
#endif
		return 0;
	}

	/*
	 * Convert back w to plain 32-bit format.
	 */
	mp_iNTT(logn, t1, igm, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = (uint32_t)mp_norm(t1[u], p);
	}

#if NTRUGEN_STATS
	stats_compute_w_success ++;
#endif
	return 1;
}

/*
 * Test whether a given polynomial is invertible modulo X^n+1 and q, for
 * a small prime q <= 1448 such that q-1 is a multiple of 256 (so q = 257
 * or q = 769). We use the NTT modulo r, which is an odd multiple of q.
 *
 * Rules:
 *    logn = 8, 9 or 10
 *    q = 257 or 769
 *    floor(x/q) = floor((x*m)/(2^j)) for all x in [0, 2^32-1]
 *    r is odd
 *    r is a multiple of q
 *    (3/2)*2^30 < r < 2^31
 *    r0i = -1/r mod 2^32
 *    s = (2^32)*t mod r, for some t such that t^128 = -1 mod q
 *
 * Returned value: 1 if invertible mod X^n+1 mod q, 0 otherwise.
 *
 * RAM usage: 3*n words in tmp[]
 */
static int
poly_is_invertible_mq(unsigned logn, const int8_t *restrict f, uint32_t q,
	uint32_t r, uint32_t r0i, uint32_t s,
	uint32_t m, unsigned j, uint32_t *restrict tmp)
{
	/*
	 * We cannot use the full-width NTT because q-1 is not a multiple
	 * of 2^(n+1). However, we can manually halve the degree a few
	 * times by considering the following:
	 * If f = f_e(X^2) + X*f_o(X^2), then let f' = f_e(X^2) - X*f_o(X^2).
	 * If f is invertible modulo q, then there exists a g such that
	 * f*g = 1. Write g = g_e(X^2) + X*g_o(X^2). Then:
	 *
	 *   f*g = (f_e*g_e + X*f_o*g_o)(X^2) + X*(f_e*g_o + f_o*g_e)(X^2)
	 *   f'*g' = (f_e*g_e + X*f_o*g_o)(X^2) - X*(f_e*g_o + f_o*g_e)(X^2)
	 *
	 * Since f*g = 1, all its odd-indexed coefficients are zero, which
	 * implies that f_e*g_o + f_o*g_e = 0. Therefore, f'*g' = f*g = 1,
	 * and f' is invertible. Thus, f*f' is invertible.
	 * Conversely if f*f' is invertible, then f and f' are invertible.
	 * We thus have an equivalent: f is invertible mod X^n+1 mod q if
	 * and only if f*f' is invertible mod X^n+1 mod q. But we can
	 * compute f*f' over the half-degree n/2, because:
	 *
	 *  f*f' = f_e^2(X^2) - (X^2)*f_o(X^2)
	 *       = (f_e^2 - X*f_o^2)(X^2)
	 *
	 * Thus, f is invertible mod X^n+1 mod q if and only if
	 * f_e^2 - X*f_o^2 is invertible mod X^(n/2)+1 mod q.
	 * We compute this degree halving transform over plain ingtegers
	 * sufficiently many times so that degree is reduced to 128, where
	 * the NTT applies.
	 *
	 * To compute over plain integers, we in fact use the usual modulus
	 * p = 2147473409, which is large enough to avoid wrap-around as
	 * long as we use reductions modulo q on non-NTT coefficients between
	 * any two halvings.
	 */
	size_t n = (size_t)1 << logn;
	uint32_t *t1 = tmp;
	uint32_t *gm = t1 + n;
	uint32_t *igm = gm + n;
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	uint32_t R2 = PRIMES[0].R2;
	mp_mkgmigm(logn, gm, igm, PRIMES[0].g, PRIMES[0].ig, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = (uint32_t)f[u];
	}
	while (logn > 7) {
		/*
		 * Perform the even/odd split and compute f_e^2 - X*f_o^2
		 * in NTT. The polynomial f_e^2 - X*f_o^2 is the Galois
		 * norm of f, and in our NTT representation it is obtained
		 * by multiplying the coefficients pairwise.
		 */
		for (size_t u = 0; u < n; u ++) {
			t1[u] = mp_set(*(int32_t *)&t1[u], p);
		}
		mp_NTT(logn, t1, gm, p, p0i);
		for (size_t u = 0; u < n; u += 2) {
			uint32_t xe = t1[u + 0];
			uint32_t xo = t1[u + 1];
			t1[u >> 1] = mp_montymul(
				mp_montymul(xe, xo, p, p0i), R2, p, p0i);
		}
		logn --;
		n >>= 1;
		mp_iNTT(logn, t1, igm, p, p0i);

		/*
		 * Convert back to plain and reduce modulo q.
		 */
		for (size_t u = 0; u < n; u ++) {
			uint32_t a = (uint32_t)mp_norm(t1[u], p);

			/*
			 * Replace a with -a - 1 if a < 0; keep sign in sa.
			 */
			uint32_t sa = tbmask(a);
			a ^= sa;

			/*
			 * Reduce a modulo q.
			 */
			a -= q * (uint32_t)(((uint64_t)a * m) >> j);

			/*
			 * Adjust the value if the source was negative.
			 */
			a += sa & (q - 1 - (a << 1));
			t1[u] = a;
		}
	}

	/*
	 * Now that degree n = 128, we can use the NTT (mod r) to test the
	 * invertibility modulo q.
	 */
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mp_set(*(int32_t *)&t1[u], r);
	}
	mp_mkgm7(gm, s, r, r0i);
	mp_NTT(logn, t1, gm, r, r0i);
	uint32_t e = 0;
	for (size_t u = 0; u < n; u ++) {
		uint32_t a = t1[u];
		a -= q * (uint32_t)(((uint64_t)a * m) >> j);
		e |= a - 1;
	}
	return 1 - (int)(e >> 31);
}

/*
 * Wrappers for q = 257 and q = 769; they are non-static for test purposes.
 */
#define poly_is_invertible_mod257   Zn(poly_is_invertible_mod257)
int
poly_is_invertible_mod257(unsigned logn,
	const int8_t *restrict f, uint32_t *restrict tmp)
{
	return poly_is_invertible_mq(logn, f, 257,
		2147483519, 2413838209, 774, 4278255361, 40, tmp);
}

#define poly_is_invertible_mod769   Zn(poly_is_invertible_mod769)
int
poly_is_invertible_mod769(unsigned logn,
	const int8_t *restrict f, uint32_t *restrict tmp)
{
	return poly_is_invertible_mq(logn, f, 769,
		2147482485, 2202878755, 16282, 2859588109, 41, tmp);
}

/*
 * Some useful NTT properties:
 * ---------------------------
 *
 * In NTT representation: given NTT(g) for g of degree n/2, then NTT(g(X^2))
 * is obtained by duplicating each coefficient:
 *    NTT(g(X^2))[2*j] = NTT(g(X^2))[2*j+1] = NTT(g)[j]
 * (This is true for our NTT coefficients, which are in bit reversal order.)
 *
 * NTT(X) is such that:
 *    NTT(X)[2*j + 0] = gm[j + n/2]
 *    NTT(X)[2*j + 1] = -gm[j + n/2]
 */

/* see ntrugen.h */
int
BAT_keygen(unsigned logn,
	int8_t *restrict f, int8_t *restrict g,
	int8_t *restrict F, int8_t *restrict G, int32_t *w,
	ntrugen_rng rng, void *restrict rng_context,
	void *restrict tmp, size_t tmp_len)
{
	/*
	 * Ensure that the tmp[] buffer has proper alignment for 64-bit
	 * access, and check the length.
	 */
	if (tmp_len < 7) {
		return -1;
	}
	if (logn < 8 || logn > 10) {
		return -1;
	}
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 7) & ~(uintptr_t)7;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint32_t *tt32 = (void *)utmp2;
	if (tmp_len < ((size_t)20 << logn)) {
		return -1;
	}

	for (;;) {
		uint32_t gamma2;
		uint32_t bound_norm2_fg;
		const ntru_profile *prof;

		/*
		 * Generate f and g. This ensures that f and g have
		 * odd parity, and are both invertible modulo q.
		 */
		switch (logn) {
		case 8:
			gauss_sample_poly(logn, f, gauss_BAT_128_256,
				rng, rng_context);
			gauss_sample_poly(logn, g, gauss_BAT_128_256,
				rng, rng_context);
			gamma2 = 1;
			bound_norm2_fg = 181;
			prof = &SOLVE_BAT_128_256;
			break;
		case 9:
			do {
				gauss_sample_poly(logn, f, gauss_BAT_257_512,
					rng, rng_context);
			} while (!poly_is_invertible_mod257(logn, f, tt32));
			do {
				gauss_sample_poly(logn, g, gauss_BAT_257_512,
					rng, rng_context);
			} while (!poly_is_invertible_mod257(logn, g, tt32));
			gamma2 = 1;
			bound_norm2_fg = 363;
			prof = &SOLVE_BAT_257_512;
			break;
		case 10:
			do {
				gauss_sample_poly(logn, f, gauss_BAT_769_1024,
					rng, rng_context);
			} while (!poly_is_invertible_mod769(logn, f, tt32));
			do {
				gauss_sample_poly(logn, g, gauss_BAT_769_1024,
					rng, rng_context);
			} while (!poly_is_invertible_mod769(logn, g, tt32));
			gamma2 = 5;
			bound_norm2_fg = 2671;
			prof = &SOLVE_BAT_769_1024;
			break;
		default:
			/*
			 * Other degrees are not supported.
			 */
			return -1;
		}

		/*
		 * Check that (g,gamma*f) has an acceptable norm. Maximum
		 * allowed norm is sqrt(n) * sigma * sqrt(gamma^2 + 1)
		 * less than 1.17*sqrt(q)).
		 */
		uint32_t norm2_fg = gamma2 * poly_sqnorm(logn, f)
			+ poly_sqnorm(logn, g);
		if (norm2_fg > bound_norm2_fg) {
			continue;
		}

		/*
		 * Solve the NTRU equation.
		 */
#if NTRUGEN_STATS
		stats_solve_attempt ++;
#endif
		int err = solve_NTRU(prof, logn, g, f, tt32);
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
		 * G and F are at the start of tt32, in that order (since
		 * we use the BAT convention here).
		 */
		size_t n = (size_t)1 << logn;
		int8_t *tG = (int8_t *)tt32;
		int8_t *tF = tG + n;
		uint32_t *tw = (uint32_t *)(tF + n);

		/*
		 * If the computation of w fails, then we have to try
		 * again (such failures are rare).
		 */
		if (!compute_w(logn, f, g, tF, tG, tw)) {
			continue;
		}

		/*
		 * Return F, G and w in the provided arrays.
		 */
		if (F != NULL) {
			memmove(F, tF, n);
		}
		if (G != NULL) {
			memmove(G, tG, n);
		}
		if (w != NULL) {
			memmove(w, tw, n * sizeof *tw);
		}

		/*
		 * Also return the values in the proper order (caller
		 * expects (F,G), not (G,F)). We also need to consider that
		 * alignment may have moved tt32 a few bytes beyond tmp.
		 */
		int8_t *dF = (int8_t *)tmp;
		int8_t *dG = dF + n;
		int32_t *dw = (int32_t *)(dG + n);
		memmove(tw + n, tG, n);
		tG = (int8_t *)(tw + n);
		memmove(dF, tF, n);
		memmove(dG, tG, n);
		if ((void *)tw != (void *)dw) {
			memmove(dw, tw, n * sizeof *tw);
		}

		return 0;
	}
}

/* see ntrugen.h */
int
BAT_recover_G(unsigned logn,
	int8_t *restrict G,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, void *tmp, size_t tmp_len)
{
	const ntru_profile *prof;
	switch (logn) {
	case 8:   prof = &SOLVE_BAT_128_256;   break;
	case 9:   prof = &SOLVE_BAT_257_512;   break;
	case 10:  prof = &SOLVE_BAT_769_1024;  break;
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

	/*
	 * We call recover_G() with -q instead of q to account for the
	 * BAT convention (f*G - g*F = -q).
	 */
	int r = recover_G(logn,
		-(int32_t)prof->q, prof->coeff_FG_limit[logn],
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

/*
 * Note: we do not provide a compute_public() function for BAT because
 * it is not needed for private key generation, and it needs a way to
 * compute inverses modulo X^n+1 and modulo q, which is cumbersome since
 * q-1 is not a multiple of 2048. A BAT implementation should already
 * have efficient code for computations modulo q, and it does not make
 * much sense to duplicate it here, in a less efficient way.
 */
