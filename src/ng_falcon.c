#include "ng_inner.h"

const ntru_profile SOLVE_Falcon_256 = {
	12289,
	2, 8,
	{ 1, 1, 2, 3, 4, 8, 14, 27, 53, 104, 207 },
	{ 1, 2, 3, 6, 11, 21, 40, 78, 155, 308 },
	{ 1, 1, 2, 2, 2, 3, 3, 4, 5, 7 },
	16,
	{ 0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};

const ntru_profile SOLVE_Falcon_512 = {
	12289,
	9, 9,
	{ 1, 1, 2, 3, 4, 8, 14, 27, 53, 104, 207 },
	{ 1, 2, 3, 6, 11, 21, 40, 78, 155, 308 },
	{ 1, 1, 2, 2, 2, 3, 3, 4, 5, 7 },
	13,
	{ 0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};

const ntru_profile SOLVE_Falcon_1024 = {
	12289,
	10, 10,
	{ 1, 1, 2, 3, 4, 8, 14, 27, 53, 104, 207 },
	{ 1, 2, 3, 6, 11, 21, 40, 78, 155, 308 },
	{ 1, 1, 2, 2, 2, 3, 3, 4, 5, 7 },
	11,
	{ 0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};

/* Falcon, q = 12289, n = 256 -> kmax = 24 */
const uint16_t gauss_Falcon_256[] = {
	24,
	    1,     3,     6,    11,    22,    40,    73,   129,
	  222,   371,   602,   950,  1460,  2183,  3179,  4509,
	 6231,  8395, 11032, 14150, 17726, 21703, 25995, 30487,
	35048, 39540, 43832, 47809, 51385, 54503, 57140, 59304,
	61026, 62356, 63352, 64075, 64585, 64933, 65164, 65313,
	65406, 65462, 65495, 65513, 65524, 65529, 65532, 65534
};

/* Falcon, q = 12289, n = 512 -> kmax = 17 */
const uint16_t gauss_Falcon_512[] = {
	17,
	    1,     4,    11,    28,    65,   146,   308,   615,
	 1164,  2083,  3535,  5692,  8706, 12669, 17574, 23285,
	29542, 35993, 42250, 47961, 52866, 56829, 59843, 62000,
	63452, 64371, 64920, 65227, 65389, 65470, 65507, 65524,
	65531, 65534
};

/* Falcon, q = 12289, n = 1024 -> kmax = 12 */
const uint16_t gauss_Falcon_1024[] = {
	12,
	    2,     8,    28,    94,   280,   742,  1761,  3753,
	 7197, 12472, 19623, 28206, 37329, 45912, 53063, 58338,
	61782, 63774, 64793, 65255, 65441, 65507, 65527, 65533
};

/* see ntrugen.h */
int
Falcon_keygen(unsigned logn,
	int8_t *restrict f, int8_t *restrict g,
	int8_t *restrict F, int8_t *restrict G,
	ntrugen_rng rng, void *restrict rng_context,
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

	const ntru_profile *prof;
	for (;;) {
		/*
		 * Generate f and g.
		 */
		switch (logn) {
		case 8:
			prof = &SOLVE_Falcon_256;
			gauss_sample_poly(logn, f, gauss_Falcon_256,
				rng, rng_context);
			gauss_sample_poly(logn, g, gauss_Falcon_256,
				rng, rng_context);
			break;
		case 9:
			prof = &SOLVE_Falcon_512;
			gauss_sample_poly(logn, f, gauss_Falcon_512,
				rng, rng_context);
			gauss_sample_poly(logn, g, gauss_Falcon_512,
				rng, rng_context);
			break;
		case 10:
			prof = &SOLVE_Falcon_1024;
			gauss_sample_poly(logn, f, gauss_Falcon_1024,
				rng, rng_context);
			gauss_sample_poly(logn, g, gauss_Falcon_1024,
				rng, rng_context);
			break;
		default:
			prof = &SOLVE_Falcon_256;
			if (logn < 2 || logn > 7) {
				return -1;
			}
			gauss_sample_poly_reduced(logn, f, gauss_Falcon_256,
				rng, rng_context);
			gauss_sample_poly_reduced(logn, g, gauss_Falcon_256,
				rng, rng_context);
			break;
		}

		/*
		 * Check that (g,-f) has an acceptable norm (norm must be
		 * less than 1.17*sqrt(q)).
		 */
		if ((poly_sqnorm(logn, f) + poly_sqnorm(logn, g)) >= 16823) {
			continue;
		}

		/*
		 * Check that f is invertible modulo X^n+1 and modulo q.
		 *
		 * We compute modulo the non-prime 2147465883, which is
		 * a multiple of 12289.
		 */
		if (!poly_is_invertible(logn, f,
			2147465883, 2763744365, 248710,
			12289, 2863078533, 45, tt32))
		{
			continue;
		}

		/*
		 * For Falcon, we need to check that the orthogonalized
		 * vector also has an acceptable norm.
		 */
		size_t n = (size_t)1 << logn;
		fxr *rt1 = (fxr *)tt32;
		fxr *rt2 = rt1 + n;
		fxr *rt3 = rt2 + n;
		vect_set(logn, rt1, f);
		vect_set(logn, rt2, g);
		vect_FFT(logn, rt1);
		vect_FFT(logn, rt2);
		vect_invnorm_fft(logn, rt3, rt1, rt2, 0);
		vect_adj_fft(logn, rt1);
		vect_adj_fft(logn, rt2);
		vect_mul_realconst(logn, rt1, fxr_of(12289));
		vect_mul_realconst(logn, rt2, fxr_of(12289));
		vect_mul_autoadj_fft(logn, rt1, rt3);
		vect_mul_autoadj_fft(logn, rt2, rt3);
		vect_iFFT(logn, rt1);
		vect_iFFT(logn, rt2);
		fxr sn = fxr_zero;
		for (size_t u = 0; u < n; u ++) {
			sn = fxr_add(sn,
				fxr_add(fxr_sqr(rt1[u]), fxr_sqr(rt2[u])));
		}
		if (!fxr_lt(sn, fxr_of_scaled32(72251709809335))) {
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
		 * Return the computed F and G.
		 */
		int8_t *tF = (int8_t *)tt32;
		int8_t *tG = tF + n;
		if (F != NULL) {
			memmove(F, tF, n);
		}
		if (G != NULL) {
			memmove(G, tG, n);
		}
		if (tt32 != tmp) {
			memmove(tmp, tt32, 2 * n);
		}

		return 0;
	}
}

#if 0
/* remove -- this should be in a Falcon implementation, with mod q code,
   not in ntrugen */

/* see ntrugen.h */
int
Falcon_recover_G(unsigned logn,
	int8_t *restrict G,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, void *tmp, size_t tmp_len)
{
	/*
	 * Ensure that the tmp[] buffer is the right size; adjust
	 * alignment.
	 */
	if (tmp_len < 1) {
		return -1;
	}
	if (logn < 2 || logn > 10) {
		return -1;
	}
	const ntru_profile *prof;
	switch (logn) {
	case 9:  prof = &SOLVE_Falcon_512;  break;
	case 10: prof = &SOLVE_Falcon_1024; break;
	default: prof = &SOLVE_Falcon_256;  break;
	}
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 1) & ~(uintptr_t)1;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint16_t *tt16 = (void *)utmp2;
	if (tmp_len < ((size_t)4 << logn)) {
		return -1;
	}

	/*
	 * f*G - g*F = q, hence G = (q + g*F)/f
	 * We work modulo q, hence G = (g*F)/f mod q.
	 * Keygen should have made sure that f is invertible modulo q.
	 */
	size_t n = (size_t)1 << logn;
	uint32_t r = 0;
	uint16_t *t1 = tt16;
	uint16_t *t2 = t1 + n;
	mq12289_poly_set_small(logn, t1, g);
	mq12289_poly_set_small(logn, t2, F);
	mq12289_NTT(logn, t1);
	mq12289_NTT(logn, t2);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mq12289_tomonty(mq12289_montymul(t1[u], t2[u]));
	}
	mq12289_poly_set_small(logn, t2, f);
	mq12289_NTT(logn, t2);
	for (size_t u = 0; u < n; u ++) {
		/* top bit of r is set if t2[u] is zero (i.e. 12289) */
		r |= 12288 - (uint32_t)t2[u];
		t2[u] = mq12289_div(t1[u], t2[u]);
	}
	mq12289_iNTT(logn, t2);
	int32_t lim = prof->coeff_FG_limit[logn];
	uint8_t *Gt = (uint8_t *)tmp;
	for (size_t u = 0; u < n; u ++) {
		int32_t z = mq12289_snorm(t2[u]);
		/* top bit of r is set if z > +lim or z < -lim */
		r |= (uint32_t)(lim - z);
		r |= (uint32_t)(lim + z);
		Gt[u] = (uint8_t)z;
	}
	if (G != NULL) {
		memmove(G, Gt, n);
	}
	if ((void *)Gt != tmp) {
		memmove(tmp, Gt, n);
	}
	return -(int)(r >> 31);
}

/* see ntrugen.h */
int
Falcon_compute_public(unsigned logn,
	uint16_t *h, const int8_t *restrict f, const int8_t *restrict g,
	void *tmp, size_t tmp_len)
{
	/*
	 * Ensure that the tmp[] buffer has proper alignment for 64-bit
	 * access.
	 */
	if (tmp_len < 1) {
		return -1;
	}
	if (logn < 2 || logn > 10) {
		return -1;
	}
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 1) & ~(uintptr_t)1;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint16_t *tt16 = (void *)utmp2;
	if (tmp_len < ((size_t)4 << logn)) {
		return -1;
	}

	size_t n = (size_t)1 << logn;
	uint16_t *t1 = tt16;
	uint16_t *t2 = t1 + n;
	mq12289_poly_set_small(logn, t1, f);
	mq12289_poly_set_small(logn, t2, g);
	mq12289_NTT(logn, t1);
	mq12289_NTT(logn, t2);
	uint32_t r = 0;
	for (size_t u = 0; u < n; u ++) {
		r |= 12288 - (uint32_t)t1[u];
		t1[u] = mq12289_div(t2[u], t1[u]);
	}
	mq12289_iNTT(logn, t1);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mq12289_unorm(t1[u]);
	}
	if (h != NULL) {
		memmove(h, t1, n * sizeof *t1);
	}
	if ((void *)t1 != tmp) {
		memmove(tmp, t1, n * sizeof *t1);
	}
	return -(int)(r >> 31);
}
#endif
