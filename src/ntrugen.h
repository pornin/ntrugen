/* ==================================================================== */
/*
 * NTRUGEN
 * -------
 *
 * This code implements the key pair generation for the BAT, Falcon and
 * Hawk cryptographic algorithms. For all three algorithm, the private
 * key mainly consists of four polynomials f, g, F, and G, with small
 * integer coefficients and considered modulo X^n+1 for a given degree n.
 * The degree n is a power of 2, and is provided as the 'logn' parameter
 * (with n = 2^logn). The polynomials are linked together through the
 * NTRU equation f*G - g*F = q, for an integer q which depends on the
 * algorithm. Each algorithm has its own unique distributions and
 * requirements for f, g, F and G. BAT private keys include a fifth
 * integer polynomial called w. The implementation follows all the
 * algorithm-specific requirements and produces f, g, F, G (and w)
 * appropriate for cryptographic usage, provided that used random source
 * is cryptographically secure.
 *
 * The code is believed to be secure and, in particular, constant-time.
 * Moreover, it uses no floating-point operations; all internal
 * computations rely only on integer operations. All the code is
 * thread-safe (there is no mutable global state).
 *
 * Each of the function defined below expects a "temporary buffer"
 * (called tmp, of size tmp_len bytes). The size of that buffer is
 * specified for each function with an expression similar to "24*n + 7";
 * internally, the code will need 8-byte alignment in memory, and the
 * extra 7 bytes are meant to ensure that such an alignment can be
 * achieved. If the caller can make sure that the tmp buffer is already
 * 8-byte aligned, then the extra 7 bytes are not needed. On some
 * architectures, enforcing a 32-byte or even 64-byte alignment _may_
 * slightly help with performance.
 *
 * The contents of the tmp buffer on input are ignored. On output, some
 * of the value produced by the function are copied in that buffer. For
 * instance, ntrugen_Falcon() expects as parameters some pointers to the
 * output buffers for f, g, F and G (n bytes each), but the F and G
 * pointers can be set to NULL because the F and G polynomials are
 * _also_ copied to the start of tmp on output, so that a caller in a
 * scarce memory situation may save a couple of kilobytes.
 *
 * All three key pair generation functions need 24*n bytes for tmp (+7 if
 * not 8-byte aligned), and an extra 2*n bytes for the f and g output
 * arrays (while F and G can be retrived from tmp on output), for a total
 * of 26624 bytes at n = 1024 (13312 bytes at n = 512).
 */

#ifndef NTRUGEN_H__
#define NTRUGEN_H__

#include <stddef.h>
#include <stdint.h>

/*
 * Type for a random source. The source is supposed to provide
 * bytes indistinguishable from uniform random selection. When this
 * function is called, it should write len random bytes into the
 * provided buffer (dst). The context pointer (ctx) is an arbitrary
 * value that is supposed to be used by the random source implementation
 * for referencing any running state that it needs.
 *
 * The random source is not supposed to fail; there is no provision for
 * returning an error code. In practice, the random source should be a
 * cryptographically secure PRNG seeded with a high-entropy value; the
 * ntrugen_prng.h file documents the API for such a PRNG.
 *
 * Optimization note: in the case of BAT and Falcon, the ntrugen code will
 * always invoke this source always by chunks of 512 bytes (i.e. len == 512).
 * The Hawk implementation makes only a few (smaller) requests.
 */
typedef void (*ntrugen_rng)(void *ctx, void *dst, size_t len);

/* ==================================================================== */
/*
 * BAT
 * ---
 *
 * BAT is a key encapsulation mechanism specified in:
 *    https://eprint.iacr.org/2022/031
 *
 * This code supports three BAT variants, selected by the degree logarithm
 * (logn parameter):
 *
 *  name            q   degree n   logn   security            tmp_len (bytes)
 *  BAT-128-256    128    256        8    80 bits               6151
 *  BAT-257-512    257    512        9    NIST-I (128 bits)    12295
 *  BAT-769-1024   769   1024       10    NIST-V (256 bits)    24583
 *
 * Private key elements are the polynomials f, g, F, G and w. All five
 * polynomials contain n elements. Elements of f and g are in [-2..+2]
 * ([-3..+3] for BAT-769-1024). Elements of F and G are in [-31..+31].
 * Elements of w are in [-65536..+65536].
 * The temporary buffer (tmp) requires 24*n bytes with 64-bit alignment;
 * providing 24*n + 7 bytes always works, regardless of the buffer alignment.
 *
 * BAT uses the convention that g*F - f*G = q; this API follows that
 * convention (equivalent, f*G - g*F = -257 for BAT-257-512).
 */

/*
 * Do key pair generation for BAT. The f, g, F, G and w polynomials are
 * filled with a new key pair.
 *    logn          degree logarithm (8, 9 or 10)
 *    f             private key polynomial (size: 2^logn)
 *    g             private key polynomial (size: 2^logn)
 *    F             private key polynomial (size: 2^logn) (may be NULL)
 *    G             private key polynomial (size: 2^logn) (may be NULL)
 *    w             private key polynomial (size: 2^logn) (may be NULL)
 *    rng           randommness source (invoked by chunks of 512 bytes)
 *    rng_context   context structure for the RNG
 *    tmp           buffer for temporary values
 *    tmp_len       length of the buffer (in bytes)
 * Returned value: 0 on success, -1 on error (invalid parameters).
 *
 * On success:
 *  - f and g are filled with n-element polynomials. The f and g parameters
 *    MUST NOT be NULL.
 *  - The corresponding F, G and w are computed and returned, in that order,
 *    in the provided tmp buffer (F and G use n bytes each, w uses n
 *    elements of type int32_t).
 *  - The F, G and w polynomials are also written into the buffers provided
 *    as parameters with those names, if not NULL. The F, G and/or w
 *    parameters may be NULL. Regardless of whether the pointers provided
 *    as parameters are NULL or not, the polynomials are always computed
 *    and copied into the tmp buffer.
 *
 * The length of the temporary buffer should be at least 24*(2^logn) + 7 bytes.
 */
int ntrugen_BAT_keygen(unsigned logn,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G, int32_t *w,
	ntrugen_rng rng, void *rng_context, void *tmp, size_t tmp_len);

/*
 * Given the private key elements f, g and F, recompute the fourth
 * polynomial G.
 *
 * The G polynomial is returned in the first n bytes of the tmp buffer.
 * If the G parameter is not NULL, then the G polynomial is also copied
 * into the buffer pointed to by G.
 *
 * Returned value: 0 on success, -1 on error (invalid parameters or key).
 * The length of the temporary buffer should be at least 12*(2^logn) + 7 bytes.
 */
int ntrugen_BAT_recover_G(unsigned logn,
	int8_t *G, const int8_t *f, const int8_t *g, const int8_t *F,
	void *tmp, size_t tmp_len);

/* ==================================================================== */
/*
 * Falcon
 * ------
 *
 * Falcon is a signature algorithm specified in:
 *    https://falcon-sign.info/falcon.pdf
 * It has been selected for standardization by NIST and should thus
 * presumably become a NIST standard at some point in the future.
 *
 * This code supports standard and reduced Falcon parameters. Falcon
 * uses q = 12289; standard parameters use degree 512 (logn = 9, for
 * Falcon-512) and 1024 (logn = 10, for Falcon-1024). This code also
 * supports reduced degrees 4 to 256 (logn = 2 to 8); these reduced
 * variants are meant for test purposes. Falcon-256 might actually provide
 * enough security for some specific use cases, especially for signature
 * with only short-term security requirements.
 *
 * Private key elements are the polynomials f, g, F and G. Elements of f
 * and g have a limited range that depends on the degree:
 *    degree n    logn    max range for (f,g)
 *    4 to 32    2 to 5   [-127..+127]
 *       64         6     [-96..+96]
 *      128         7     [-48..+48]
 *      256         8     [-24..+24]
 *      512         9     [-17..+17]
 *     1024        10     [-12..+12]
 * Elements of F and G all fit in [-127..+127].
 * The temporary buffer (tmp) requires 24*n bytes with 64-bit alignment;
 * providing 24*n + 7 bytes always works, regardless of the buffer alignment.
 */

/*
 * Do key pair generation for Falcon. The f, g, F and G polynomials are
 * filled with a new key pair.
 *    logn          degree logarithm (8, 9 or 10)
 *    f             private key polynomial (size: 2^logn)
 *    g             private key polynomial (size: 2^logn)
 *    F             private key polynomial (size: 2^logn) (may be NULL)
 *    G             private key polynomial (size: 2^logn) (may be NULL)
 *    rng           randommness source (invoked by chunks of 512 bytes)
 *    rng_context   context structure for the RNG
 *    tmp           buffer for temporary values
 *    tmp_len       length of the buffer (in bytes)
 * Returned value: 0 on success, -1 on error (invalid parameters).
 *
 * On success:
 *  - f and g are filled with n-element polynomials. The f and g parameters
 *    MUST NOT be NULL.
 *  - The corresponding F and G are computed and returned, in that order,
 *    in the provided tmp buffer (F and G use n bytes each).
 *  - The F and G polynomials are also written into the buffers provided
 *    as parameters with those names, if not NULL. The F and/or G
 *    parameters may be NULL. Regardless of whether the pointers provided
 *    as parameters are NULL or not, the polynomials are always computed
 *    and copied into the tmp buffer.
 *
 * The length of the temporary buffer should be at least 24*(2^logn) + 7 bytes.
 */
int ntrugen_Falcon_keygen(unsigned logn,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G,
	ntrugen_rng rng, void *rng_context, void *tmp, size_t tmp_len);

/* ==================================================================== */
/*
 * Hawk
 * ----
 *
 * Hawk is a signature algorithm specified in:
 *    https://eprint.iacr.org/2022/1155
 *
 * This code supports three sets of parameters for Hawk:
 *    name        degree n   logn   (f,g) range   security
 *    Hawk-256       256       8      [-5..+5]    "challenge" (~64 bits)
 *    Hawk-512       512       9      [-6..+6]    NIST-I (128 bits)
 *    Hawk-1024     1024      10      [-9..+9]    NIST-V (256 bits)
 *
 * Private key elements are the polynomials f, g, F and G. Elements of f
 * and g have a limited range that depends on the degree (see table
 * above). Elements of F and G all fit in [-127..+127].
 * The temporary buffer (tmp) requires 24*n bytes with 64-bit alignment;
 * providing 24*n + 7 bytes always works, regardless of the buffer alignment.
 */

/*
 * Do key pair generation for Hawk. The f, g, F, G, q00, q01 and q11
 * polynomials are filled with a new key pair.
 *    logn          degree logarithm (8, 9 or 10)
 *    f             private key polynomial (size: 2^logn)
 *    g             private key polynomial (size: 2^logn)
 *    F             private key polynomial (size: 2^logn) (may be NULL)
 *    G             private key polynomial (size: 2^logn) (may be NULL)
 *    q00           public key polynomial (size: 2^logn) (may be NULL)
 *    q01           public key polynomial (size: 2^logn) (may be NULL)
 *    q11           public key polynomial (size: 2^logn) (may be NULL)
 *    seed          seed for (f,g) re-gen (size: 8 + 2^(logn-5)) (may be NULL)
 *    rng           randommness source (invoked by chunks of 512 bytes)
 *    rng_context   context structure for the RNG
 *    tmp           buffer for temporary values
 *    tmp_len       length of the buffer (in bytes)
 * Returned value: 0 on success, -1 on error (invalid parameters).
 *
 * On success:
 *  - f and g are filled with n-element polynomials. The f and g parameters
 *    MUST NOT be NULL.
 *  - The corresponding F, G, q00, q01, q11 and seed are computed and returned,
 *    in that order, in the provided tmp buffer (F and G have type int8_t[n],
 *    q00 and q01 have type int16_t[n], q11 has type int32_t[n], the seed
 *    is uint8_t[seed_len]).
 *  - The F, G, q00, q01 and q11 polynomials are also written into the buffers
 *    provided as parameters with those names, if not NULL. The F, G, q00,
 *    q01 and/or q11 parameters may be NULL. Regardless of whether the
 *    pointers provided as parameters are NULL or not, the polynomials are
 *    always computed and written into the tmp buffer.
 *  - The seed used to generate (f,g) is written into the provided seed
 *    buffer if it is not NULL. The seed length is 8 + 2^(logn-5) bytes
 *    (i.e. 16, 24 or 40 bytes, for logn = 8, 9 or 10). The seed is also
 *    always written to tmp[] right after q11. The seed can be used to
 *    regenerate the same (f,g) with ntrugen_Hawk_regen_fg().
 *
 * The length of the temporary buffer should be at least 24*(2^logn) + 7 bytes.
 */
int ntrugen_Hawk_keygen(unsigned logn,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G,
	int16_t *q00, int16_t *q01, int32_t *q11,
	void *seed, ntrugen_rng rng, void *rng_context,
	void *tmp, size_t tmp_len);

/*
 * Regenerate (f,g) from the given seed. The seed length is 16, 24 or 40
 * bytes, for logn = 8, 9 or 10. Note that no validations are performed
 * on the returned (f,g).
 */
void ntrugen_Hawk_regen_fg(unsigned logn,
	int8_t *f, int8_t *g, const void *seed);

/*
 * Given the private key elements f, g and F, recompute the fourth
 * private polynomial G.
 *
 * The G polynomial is returned in the first n bytes of the tmp buffer.
 * If the G parameter is not NULL, then the G polynomial is also copied
 * into the buffer pointed to by G.
 *
 * Returned value: 0 on success, -1 on error (invalid parameters or key).
 * The length of the temporary buffer should be at least 12*(2^logn) + 7 bytes.
 */
int ntrugen_Hawk_recover_G(unsigned logn,
	int8_t *G, const int8_t *f, const int8_t *g, const int8_t *F,
	void *tmp, size_t tmp_len);

/*
 * Given the private key elements f, g, F and G, recompute the public
 * polynomials q00, q01 and q11.
 *
 * The recomputed polynomials are stored into the provided tmp buffer,
 * in that order:
 *    q00   int16_t[n]
 *    q01   int16_t[n]
 *    q11   int32_t[n]
 * If the q00, q01 and/or q11 parameters are not NULL, then the corresponding
 * polynomials are also copied into these buffers.
 *
 * Returned value: 0 on success, -1 on error (invalid parameters or key).
 * The length of the temporary buffer should be at least 20*(2^logn) + 7 bytes.
 */
int ntrugen_Hawk_recover_qq(unsigned logn,
	int16_t *q00, int16_t *q01, int32_t *q11,
	const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G,
	void *tmp, size_t tmp_len);

#endif
