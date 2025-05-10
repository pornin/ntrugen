# This is a Python implementation of the key pair generation of Falcon.
# THIS CODE IS FOR TEST/DEBUG PURPOSES ONLY. It is meant to show the
# computations performed by the ntrugen C code, in a format more amenable
# to specification and printing out intermediate values. It may leak
# private information through side channels; in particular, it is absolutely
# _not_ constant-time.
#
# This implementation should follow the same steps as the C code, and
# produce the exact same key pair for the source seed.
#
# Example usage (for degree n = 512, i.e. logn = 9):
#
#    prng = SHAKE256x4(b'test')
#    (f, g, F, G, h) = Falcon_keygen(9, prng)
#
# Equivalent C code:
#
#    ntrugen_prng_shake_context rng;
#    const char *seed = "test";
#    uint8_t tmp[(24 << 9) + 7];
#    int8_t f[512], g[512], F[512], G[512];
#
#    ntrugen_prng_shake_init(&rng, seed, strlen(seed));
#    if (ntrugen_Falcon_keygen(9, f, g, F, G,
#        ntrugen_prng_shake_out, &rng, tmp, sizeof tmp) < 0)
#    {
#        abort();
#    }
#
# (The C code does not compute the public key h = g/f mod x^n+1 mod q,
# because the corresponding code is more naturally part of the verification
# implementation, which is not provided in ntrugen. The Python code below
# includes the computation of the public key.)

# We use SHA-256 for consistency check of generated tables of constants,
# and SHAKE256 for the source PRNG.
import hashlib
import math

# Target small prime for Falcon.
q = 12289

# We need a primality test to generate the constant small primes p_i;
# all values on which the test will run will be less than 2^31. We
# use Erathosthenes's sieve to get all small primes less than sqrt(2^31),
# then each primality test uses trial divisions.

# primetest_pp[] receives all prime integers lower than 46341.
primetest_pp = []

# Initialization function: compute the contents of primetest_pp[].
def mk_primetest():
    max_p = 46431   # ceil(sqrt(2^31))
    max_q = 216     # ceil(sqrt(max_p))
    tab = [True] * max_p
    tab[0] = False
    tab[1] = False
    for i in range(2, max_q):
        if tab[i]:
            for j in range(2*i, max_p, i):
                tab[j] = False
    for i in range(0, max_p):
        if tab[i]:
            primetest_pp.append(i)

mk_primetest()

# Check whether the provided integer is prime. The parameter MUST be
# non-negative and lower than 2^31.
def is_small_prime(p):
    assert p >= 0 and p < 2**31
    if p <= 2:
        return p == 2
    if (p & 1) == 0:
        return False
    for j in range(0, len(primetest_pp)):
        t = primetest_pp[j]
        if t * t > p:
            return True
        if p % t == 0:
            return False
    return True

# Let (p_i) be the list of small primes p such that p < 2^31 and
# p = 1 mod 2048. We number them in reverse order; thus, p_0 = 2147473409,
# p_1 = 2147389441 and p_2 = 2147387393.
# For each depth d (0 <= d <= 10), we define the "small" modulus and the
# "large" modulus as products of p_i for i = 0 up to some limit:
#   mod_small[j] = \prod_{i=0}^{mod_small_bl[j]-1} p_i
#   mod_large[j] = \prod_{i=0}^{mod_large_bl[j]-1} p_i
#   mod_small_inc[j] = mod_small[j] * p_{mod_small[j]}
mod_small_bl = [ 1, 1, 2, 3, 4, 8, 14, 27, 53, 104, 207 ]
mod_large_bl = [ 1, 2, 3, 6, 11, 21, 40, 78, 155, 308 ]
mod_small = []
mod_large = []
mod_small_inc = []

def mk_moduli():
    # We hash the concatenation of all used p_i to check against a
    # reference value for consistency checks.
    h = hashlib.sha256()

    # mod_small[j] = \prod_{i < mod_small[j]} p_i
    # mod_large[j] = \prod_{i < mod_large[j]} p_i
    p = 2**31 + 1
    m = 1
    i = 0
    while i < mod_large_bl[len(mod_large_bl) - 1]:
        p -= 2048
        if is_small_prime(p):
            m *= p
            i += 1
            while len(mod_small) < len(mod_small_bl) and mod_small_bl[len(mod_small)] == i:
                mod_small.append(m)
            while len(mod_small_inc) < len(mod_small_bl) and mod_small_bl[len(mod_small_inc)] + 1 == i:
                mod_small_inc.append(m)
            while len(mod_large) < len(mod_large_bl) and mod_large_bl[len(mod_large)] == i:
                mod_large.append(m)
            h.update(int(p).to_bytes(4, byteorder='little'))
    assert len(mod_small) == 11
    assert len(mod_small_inc) == 11
    assert len(mod_large) == 10
    assert h.hexdigest() == '383babc06a50dafbc446cecc9043ee30ac56d598e3ef057f4ab88e0ef0368f54'

mk_moduli()
p0 = 2147473409
p1 = 2147389441

# word_win[] contains the number of top 31-bit words we have to consider on
# coefficients when doing Babai's reduction.
word_win = [ 1, 1, 2, 2, 2, 3, 3, 4, 5, 7 ]

# We use polynomials with degree less than n, taken modulo x^n+1, for
# n a power of 2 between 1 and 1024. Coefficients are integers. We
# represent a polynomial f with an array of exactly n integers, such
# that f = \sum_{i=0}^{n-1} f[i]*x^i
#
# All computations with such polynomials actually consider the
# coefficients to be within a given limited range; operations on
# coefficients are really performed modulo a given big integer m. The
# modulus m to use depends on the involved polynomial and its degree;
# possible moduli are the precomputed values in mod_small[], mod_large[]
# and mod_small_inc[], as well as some powers of two of the form
# 2^(31*k). An integer modulo m is normalized to a plain, signed integer
# by using its (unique) representant in the [-m/2..+m/2[ range. All the
# algorithm is based on the heuristic assumption that the coefficients
# are almost always lower than m/2 (in absolute value), and thus the
# exact integer value is obtained; when that assumption is not meant,
# the algorithm ultimately fails, and a new (f,g) pair is generated.
#
# This implementation is not prescriptive on how the computations are
# performed. The use of the specific defined moduli m is intended to
# help with a fast and constant-time implementation strategy in which
# coefficients are mostly in RNS representation (for integer z modulo m,
# the RNS representation is the list of (z mod p_i), for all p_i that
# divide m); then, polynomial operations can be performed with
# coefficients in the finite field GF(p_i). Each p_i was selected such
# that the NTT could be applied to the polynomials modulo x^n+1 and
# modulo p_i. For some operations, values are converted back to "plain
# integers" which are them assumed to use a given, fixed number of
# 31-bit "limbs" (with two's complement for negative value), hence the
# use of moduli 2^(31*k).
#
# In this example implementation, we merely use "plain integers" as
# provided by the Python runtime; this is inefficient (and not
# constant-time), but it makes for a simpler description of the
# algorithm, and allows easy printing of intermediate values. The
# reduction modulo m is applied at the appropriate places, in order to
# always compute exactly the same values as the reference C code, which
# uses the RNS+NTT strategy.

# Generate a new key pair, with the provided degree n. Parameter logn is
# such that n = 2^logn; it must hold that 2 <= logn <= 10 (standard
# Falcon uses logn = 9 or 10; lower values are meant for easier debug and
# tests). The 'prng' input is an object that can generate (pseudo)random
# blocks of 512 bytes, by invoking its next_512() method. The algorithm
# is fully deterministic, for a given degree n and stream of input 512-byte
# blocks. The number of requested blocks can vary.
#
# Output is (f, g, F, G, h), with h being the public key.
def Falcon_keygen(logn, prng):
    assert 2 <= logn and logn <= 10
    n = 1 << logn

    # Key pair generation consists in:
    #   1. Generate a new (f, g) pair with a Gaussian distribution.
    #   2. If the (f, g) pair is not appropriate, loop to 1. An appropriate
    #      pair must fulfill a few conditions, e.g. its L2 norm must not
    #      be too high.
    #   3. Try to find some short polynomials (F, G) such that f*G - g*F = q
    #      (modulo x^n+1). If no such (F, G) pair is found, loop to 1.
    #   4. The private key is (f, g, F, G). The public key is
    #      h = g/f mod x^n+1 mod q.
    # This function returns the private key only; computation of the public
    # key would most naturally be done with the NTT modulo q, using code
    # which is part of the Falcon signature verification (not included here).
    # The key generation process still makes sure that the public key is
    # computable (i.e. that f is invertible modulo x^n+1 and modulo q).
    while True:
        # Generate candidate (f, g) (Gaussian distribution, both with
        # odd parity).
        f = gauss_sample(logn, prng)
        g = gauss_sample(logn, prng)

        # Ensure that ||(g, -f)|| < 1.17*sqrt(q). We compute the squared
        # norm; (1.17*sqrt(q))^2 = 16822.4121
        assert math.ceil((1.17**2)*q) == 16823
        sn = 0
        for i in range(0, n):
            sn += f[i]*f[i] + g[i]*g[i]
        if sn >= 16823:
            continue

        # f must be invertible modulo x^n + 1 and modulo q. This can be
        # checked by extracting the NTT of f modulo q; f is invertible
        # if and only if none of the NTT coefficients is 0.
        qf = []
        for i in range(0, n):
            qf.append(f[i] % q)
        mq_NTT(logn, qf)
        invertible = True
        for i in range(0, n):
            if qf[i] == 0:
                invertible = False
                break
        if not(invertible):
            continue
        # We keep 1/f in NTT representation, it will be handy to compute
        # the public key later on.

        # For (f, g) to be acceptable in a Falcon key pair, the
        # orthogonalized vector must also have an acceptable norm.
        fx = []
        gx = []
        for i in range(0, n):
            fx.append(fxr_from_int(f[i]))
            gx.append(fxr_from_int(g[i]))
        vect_FFT(logn, fx)
        vect_FFT(logn, gx)
        t3 = vect_invnorm_fft(logn, fx, gx, 0)
        fx = vect_adj_fft(logn, fx)
        gx = vect_adj_fft(logn, gx)
        fx = vect_mul_realconst(logn, fx, fxr_from_int(q))
        gx = vect_mul_realconst(logn, gx, fxr_from_int(q))
        fx = vect_mul_selfadj_fft(logn, fx, t3)
        gx = vect_mul_selfadj_fft(logn, gx, t3)
        vect_iFFT(logn, fx)
        vect_iFFT(logn, gx)
        sn = fxr_from_int(0)
        for i in range(0, n):
            sn = fxr_add(sn, fxr_add(fxr_sqr(fx[i]), fxr_sqr(gx[i])))
        if sn >= fxr_from_scaled32(72251709809335):
            continue

        # (f, g) passed the basic tests; we now attempt the expensive
        # part, which is solving the NTRU equation.
        r = solve_NTRU(logn, f, g)
        if r is False:
            continue
        (F, G) = r

        # We generated the key.
        # We compute the public key h = g/f mod x^n + 1 mod q;
        # we already have f in NTT representation (in qf[]).
        qg = []
        for i in range(0, n):
            qg.append(g[i] % q)
        mq_NTT(logn, qg)
        h = []
        for i in range(0, n):
            h.append(mq_div(qg[i], qf[i]))
        mq_iNTT(logn, h)

        return (f, g, F, G, h)

# Look-up tables for the Gaussian distribution of coefficients of
# (f, g), for degrees 256, 512 and 1024. The table for degree 256
# is also used for lower degrees (with multiple samples added together).
gauss_tab8 = [
            1,     3,     6,    11,    22,    40,    73,   129,
          222,   371,   602,   950,  1460,  2183,  3179,  4509,
         6231,  8395, 11032, 14150, 17726, 21703, 25995, 30487,
        35048, 39540, 43832, 47809, 51385, 54503, 57140, 59304,
        61026, 62356, 63352, 64075, 64585, 64933, 65164, 65313,
        65406, 65462, 65495, 65513, 65524, 65529, 65532, 65534
]
gauss_tab9 = [
            1,     4,    11,    28,    65,   146,   308,   615,
         1164,  2083,  3535,  5692,  8706, 12669, 17574, 23285,
        29542, 35993, 42250, 47961, 52866, 56829, 59843, 62000,
        63452, 64371, 64920, 65227, 65389, 65470, 65507, 65524,
        65531, 65534
]
gauss_tab10 = [
            2,     8,    28,    94,   280,   742,  1761,  3753,
         7197, 12472, 19623, 28206, 37329, 45912, 53063, 58338,
        61782, 63774, 64793, 65255, 65441, 65507, 65527, 65533
]

# Generate a candidate (f, g) pair from a given pseudorandom generator.
# The source PRNG is invoked to produce only chunks of 512 bytes; each
# such chunk is then interpreted as 256 successive 16-bit unsigned integers
# (little-endian convention for each 16-bit integer), which are then used
# to produce small integers with a Gaussian distribution. The distribution
# depends on the target degree.
#
# If the generated polynomial does not have odd parity (i.e. the sum of its
# coefficients is an even intger), then the process loops and generates
# another polynomial, and so on, until an odd-parity polynomial is obtained.
# The amount of data obtained from the PRNG thus varies, but is always a
# multiple of 512 bytes.
def gauss_sample(logn, prng):
    # We use a lookup table that depends on the target degree. For
    # degree less than 8, we use the table for degree-8 but also
    # perform multiple lookups per coefficients, that we add together.
    #
    # tab = table to use for lookups
    # g = number of sampled values to add together for each output coefficient
    if logn == 10:
        tab = gauss_tab10
        g = 1
    elif logn == 9:
        tab = gauss_tab9
        g = 1
    else:
        assert logn >= 2 and logn <= 8
        tab = gauss_tab8
        g = 1 << (8 - logn)

    ru = []
    ru_ptr = 0
    while True:
        # Perform all lookups and accumulate results.
        f = []
        parity = 0
        for i in range(0, 1 << logn):
            while True:
                v = 0
                for j in range(0, g):
                    # Get next pseudorandom 16-bit integer. If the ru[]
                    # buffer is consumed, get a new one.
                    if ru_ptr >= len(ru):
                        ru = []
                        bb = prng.next_512()
                        for i in range(0, 256):
                            ru.append(bb[2 * i] + (bb[2 * i + 1] << 8))
                        ru_ptr = 0
                    y = ru[ru_ptr]
                    ru_ptr += 1

                    # Use the table: the table has 2*kmax elements; we
                    # start at -kmax, and add 1 for element table element
                    # which is strictly smaller than the random source value.
                    v -= len(tab) >> 1
                    for k in range(0, len(tab)):
                        v += ((tab[k] - y) >> 31) & 1

                # At very low degrees the value can exceed 127 in absolute
                # value; this cannot happen for the "normal" degrees
                # (256, 512 or 1024).
                if v < -127 or v > +127:
                    continue
                f.append(v)
                parity = (parity + v) & 1
                break

        # If parity is odd, we can return.
        if parity != 0:
            return f

# Given degree n (with n = 2^logn), integer polynomials f and g of
# degree less than n, and integer q, find (F, G) such that f*G - g*F = q.
# On success, F and G are returned. The process may fail for a variety of
# reasons, in which case False is returned.
def solve_NTRU(logn, f, g):
    n = 1 << logn

    # Perform the recursive solving.
    r = solve_NTRU_rec(logn, f, g, 0)
    if r is False:
        return False

    # If the recursive process found a solution, check that it fits
    # in the target range for coefficients. For Falcon, coefficients of
    # F and G must be at most 127 (in absolute value).
    (F, G) = r
    for i in range(0, n):
        if F[i] < -127 or F[i] > 127 or G[i] < -127 or G[i] > 127:
            return False

    # At depth 0, solve_NTRU_rec() checked that the computed (F, G)
    # match the NTRU equation modulo p_0. Since all coefficients of
    # f, g, F and G are at most 127 in absolute value, any
    # coefficient of f*G - g*F mod x^n+1 can be at most
    # 2*1024*127^2 = 33032192, which is much lower than p_0/2.
    # Therefore, the test modulo p_0 is sufficient to guarantee that
    # the returned solution is correct.
    return (F, G)

# Recursive NTRU solving. Provided polynomials have degree less than
# n = 2^(logn_top - depth). Output is (F, G), or False on failure.
def solve_NTRU_rec(logn_top, f, g, depth):
    # At the deepest level, f and g have degree 0, i.e. they really
    # are plain integers.
    if logn_top == depth:
        # Modulus to use.
        ms = mod_small[depth]
        hms = ms >> 1

        # We extract f and g as integers. The caller normalized them
        # modulo ms as signed integers, but since we know they should be
        # non-negative we renormalize them as unsigned.
        f = f[0] % ms
        g = g[0] % ms

        # Normally f and g are both odd, since we started with polynomials
        # with odd parity. We might get some even value if our integers
        # overflowed, in which case we report an error.
        if (f & 1) == 0 or (g & 1) == 0:
            return False

        # Get Bézout coefficients u and v such that f*u - g*v = 1.
        # We need 0 < u <= g and 0 < v < f. If f and g are not
        # co-prime then we return False (failure). This means that
        # u = 1/f mod g and v = -1/g mod f (the solution, if it exists,
        # is unique).
        #
        # We use here a binary GCD variant, with 6 variables with the
        # following initial values:
        #    a = f    u0 = 1    v0 = 0
        #    b = g    u1 = g    v1 = f - 1
        # Invariants:
        #    a = f*u0 - g*v0
        #    b = f*u1 - g*v1
        #    0 <= a <= f
        #    1 <= b <= g
        #    0 <= u0 <= g
        #    0 <= v0 < f
        #    0 <= u1 <= g
        #    0 <= v1 < f
        #
        # At each iteration:
        #    if a is odd:
        #        if a < b then swap a and b
        #        subtract b from a
        #    divide a by 2
        # We also apply the corresponding operations on u0, v0, u1 and v1 in
        # order to maintain the invariants; u0 and u1 are computed modulo g,
        # and v0 and v1 are computed modulo f. At every iteration, b is odd.
        # Algorithm stops when a reaches 0, at which point b contains the
        # GCD of f and g, and (u1,v1) are the solution (if that GCD is 1).
        a = f
        u0 = 1
        v0 = 0
        b = g
        u1 = g
        v1 = f - 1
        while a != 0:
            if (a & 1) != 0:
                if a < b:
                    (a, u0, v0, b, u1, v1) = (b, u1, v1, a, u0, v0)
                a -= b
                u0 -= u1
                v0 -= v1
                if u0 < 0:
                    u0 += g
                if v0 < 0:
                    v0 += f
            a >>= 1
            if (u0 & 1) != 0:
                u0 += g
            if (v0 & 1) != 0:
                v0 += f
            u0 >>= 1
            v0 >>= 1
        (u, v) = (u1, v1)

        assert 0 <= u and u <= g
        assert 0 < v and v < f
        assert f*u - g*v == b
        if b != 1:
            return False

        # Return (F, G) = (v*q, u*q), but return False (failure) if either
        # product is off-range. F and G should be polynomials.
        v *= q
        u *= q
        if u >= 1 << (31 * mod_small_bl[depth]) or v >= 1 << (31 * mod_small_bl[depth]):
            return False
        (F, G) = ([v], [u])
        return (F, G)

    # We are not at the deepest level. We apply the field norm and use
    # a recursive call.

    # Input polynomials are modulo x^n + 1. Recursive call will work with
    # polynomials modulo x^hn + 1, with hn = n/2.
    logn = logn_top - depth
    n = 1 << logn
    hn = n >> 1

    # f' = N(f) = fe^2 - fo^2  (mod x^(n/2)+1)
    # g' = N(g) = ge^2 - go^2  (mod x^(n/2)+1)
    # Coefficients of f' and g' use the modulus from the deeper level.
    msp = mod_small[depth + 1]
    fp = poly_field_norm(logn, f, msp)
    gp = poly_field_norm(logn, g, msp)
    r = solve_NTRU_rec(logn_top, fp, gp, depth + 1)
    if r is False:
        return False
    (Fp, Gp) = r

    # We have F' and G' from the deeper level; we then compute our own
    # (F, G) solution as:
    #    F = fx * G'(x^2)
    #    G = gx * F'(x^2)
    # with:
    #    fx = fe(x^2) - x*fo(x^2)
    #    gx = ge(x^2) - x*go(x^2)
    # The obtained F and G are unreduced at this point, and their
    # coefficients are computed modulo the "large" modulus.
    slen = mod_small_bl[depth]
    llen = mod_large_bl[depth]
    ms = mod_small[depth]
    ml = mod_large[depth]
    fs = poly_conj(logn, f, ms)
    gs = poly_conj(logn, g, ms)
    F = poly_mul(logn, gs, poly_expand(logn - 1, Fp, ml), ml)
    G = poly_mul(logn, fs, poly_expand(logn - 1, Gp, ml), ml)

    # Reduce (F, G) with Babai's round-off. This part uses fixed-point
    # approximations and computations with 64-bit integers. The steps
    # should be followed exactly if exact reproducibility of key pairs
    # is required, and they should match the C implementation.

    # Normal case is for depth >= 1; the case depth = 0 is handled specially.
    if depth >= 1:
        # Extract approximations of f and g as fixed-point values. Some
        # scaling is applied to ensure that the computation of
        # f*adj(f) + g*adj(g) (on the approximations) does not overflow;
        # the downscaling (expressed in bits) is the sum of scale_fg (a
        # public value) and scale_x (a value that depends on the actual
        # coefficients of f and g, and is thus secret).
        rlen = word_win[depth]
        blen = max(0, slen - rlen)
        scale_fg = 31*blen
        scale_FG = 31*llen
        scale_xf = poly_max_bitlength(logn, f, blen, rlen)
        scale_xg = poly_max_bitlength(logn, g, blen, rlen)
        scale_x = max(scale_xf, scale_xg)
        scale_t = min(scale_x, 15 - logn)
        scdiff = scale_x - scale_t
        fx = poly_to_fxr(logn, f, blen, rlen, scdiff)
        gx = poly_to_fxr(logn, g, blen, rlen, scdiff)

        # fx <- adj(f)/(f*adj(f) + g*adj(g))
        # gx <- adj(g)/(f*adj(f) + g*adj(g))
        vect_FFT(logn, fx)
        vect_FFT(logn, gx)
        tx = vect_norm_fft(logn, fx, gx)
        fx = vect_mul2e(logn, fx, scale_t)
        gx = vect_mul2e(logn, gx, scale_t)
        for i in range(0, hn):
            fx[i] = fxr_div(fx[i], tx[i])
            fx[i + hn] = fxr_div(fxr_neg(fx[i + hn]), tx[i])
            gx[i] = fxr_div(gx[i], tx[i])
            gx[i + hn] = fxr_div(fxr_neg(gx[i + hn]), tx[i])

        # We can now reduce (F, G) by repeatedly applying a loop that
        # "skims off" the top bits. Each iteration is assumed to gain
        # an amount of bits that depends on the top degree (a lower value
        # decreases the rejection rate but increases the number of
        # iterations; values below where measured to be a good trade-off).
        if logn_top == 10:
            reduce_bits = 11
        elif logn_top == 9:
            reduce_bits = 13
        else:
            reduce_bits = 16

        # For any iteration, coefficients of (F, G) are plain integers that
        # fit on llen limbs. This is true at this point since the unreduced
        # (F, G) were computed modulo ml, which is the product of llen small
        # primes p_i; each p_i is lower than 2^31, and thus ml < 2^(31*llen).
        FGlen = llen
        while True:
            # Extract approximations of (F, G) with the right scaling.
            tlen = scale_FG // 31
            toff = scale_FG % 31
            Fx = poly_to_fxr(logn, F, tlen, FGlen - tlen, scale_x + toff)
            Gx = poly_to_fxr(logn, G, tlen, FGlen - tlen, scale_x + toff)

            # k <- round((F*adj(f) + G*adj(g))/(f*adj(f) + g*adj(g)))
            vect_FFT(logn, Fx)
            vect_FFT(logn, Gx)
            t1 = vect_mul_fft(logn, Fx, fx)
            t2 = vect_mul_fft(logn, Gx, gx)
            t3 = vect_add(logn, t1, t2)
            vect_iFFT(logn, t3)
            k = []
            for i in range(0, n):
                k.append(fxr_round32(t3[i]))
            # Note that coefficients of k[] are in [-1518500250..+1518500250]
            # (see analysis in comments of vect_iFFT) regardless of previous
            # computations.

            # When extracting approximations, we applied scaling:
            #   (f, g) were scaled by scale_fg + scale_x
            #   (F, G) were scaled by scale_FG + scale_x
            # Thus, k is scaled by scale_FG - scale_fg, which is a public
            # value.
            # We subtract (k*f, k*g) from (F, G), taking scaling into account.
            #
            # In order to follow the exact same steps as the C code, we
            # apply the following rules:
            #
            #  - If depth > 1 and logn >= 4, then (k*f, k*g) are computed
            #    modulo mod_small_inc[depth]; the result is then subtracted
            #    (with scaling) from (F, G) modulo 2^(31*FGlen) (the C
            #    code uses RNS+NTT to compute (k*f, k*g), then converts
            #    back to plain integers, with which the shift and subtraction
            #    is performed).
            #
            #  - If depth > 1 and logn < 4, then (k*f, k*g) are computed and
            #    subtracted from (F, G) coefficient by coefficient, modulo
            #    2^(31*FGlen) (the C code uses the basic quadratic algorithm,
            #    and works with plain integers).
            #
            #  - If depth == 1, then both the multiplication (k*f, k*g) and the
            #    subtraction from (F, G) are done modulo either p_0
            #    (if FGlen = 1) or p_0*p_1 (if FGlen = 2) (the C code does
            #    the whole computation in RNS+NTT).
            #
            #  - The case depth == 0 is handled specially elsewhere. In this
            #    part of the code, we necessarily have depth == 1.
            if depth > 1:
                m2 = 1 << (31*FGlen)
                if logn >= 4:
                    m1 = mod_small_inc[depth]
                else:
                    m1 = m2
            elif depth == 1:
                if FGlen == 1:
                    m1 = p0
                else:
                    m1 = p0 * p1
                m2 = m1
            else:
                raise Exception("impossible")

            scale_k = scale_FG - scale_fg
            fk = poly_mul(logn, f, k, m1)
            gk = poly_mul(logn, g, k, m1)
            for i in range(0, n):
                F[i] = mred(F[i] - (fk[i] << scale_k), m2)
                G[i] = mred(G[i] - (gk[i] << scale_k), m2)

            # Adjust the expected length of (F, G) coefficients. Note that we
            # do _not_ test the actual lengths; for constant-time processing,
            # we only suppose that (F, G) were actually reduced. If this is
            # wrong, then this will almost surely lead to a reduction failure.
            # The final test on the NTRU equation reliably detects failures.
            if scale_FG <= scale_fg:
                break
            scale_FG = max(scale_fg, scale_FG - reduce_bits)
            while 31*(FGlen - slen) > (scale_FG - scale_fg + 30):
                FGlen -= 1
                m = 1 << (31*FGlen)
                for i in range(0, n):
                    F[i] = mred(F[i], m)
                    G[i] = mred(G[i], m)

    else:
        # Special treatment for depth == 0. This is done separately because
        # the C code does it with specialized code, making more operations
        # with integers in RNS, for performance reasons (faster and uses
        # less RAM).

        # At depth 0, slen == llen == 1, ml == ms == p0.
        # We compute F*adj(f) + G*adj(g) and f*adj(f) + g*adj(g) over
        # integers modulo p0.
        fa = poly_adj(logn, f, p0)
        ga = poly_adj(logn, g, p0)
        t1 = poly_mul(logn, F, fa, p0)
        t2 = poly_mul(logn, G, ga, p0)
        t3 = poly_add(logn, t1, t2, p0)
        t4 = poly_mul(logn, f, fa, p0)
        t5 = poly_mul(logn, g, ga, p0)
        t6 = poly_add(logn, t4, t5, p0)

        # To divide t3 by t6, we first convert them to fixed-point with
        # some down-scaling to avoid overflows.
        tx3 = []
        tx6 = []
        for i in range(0, n):
            tx3.append(fxr_from_scaled32(t3[i] << 22))
            tx6.append(fxr_from_scaled32(t6[i] << 22))
        vect_FFT(logn, tx3)
        vect_FFT(logn, tx6)
        kx = vect_div_selfadj_fft(logn, tx3, tx6)
        vect_iFFT(logn, kx)
        k = []
        for i in range(0, n):
            k.append(fxr_round32(kx[i]))

        # Compute (F - k*f, G - k*g) modulo p0.
        F = poly_sub(logn, F, poly_mul(logn, f, k, p0), p0)
        G = poly_sub(logn, G, poly_mul(logn, g, k, p0), p0)

    # We check that the solution so far seems correct, by checking that
    # f*G - g*F = q  mod x^n+1 mod p_0 (with p_0 = 2147473409). This
    # check is skipped at depth 1 because the C code omits it at that
    # depth (in the C code, that test would be inconvenient to do for
    # some memory allocation reasons); the check will still be done at
    # depth 0.
    if depth != 1:
        (tf, tg, tF, tG) = ([], [], [], [])
        for i in range(0, n):
            tf.append(f[i] % p0)
            tg.append(g[i] % p0)
            tF.append(F[i] % p0)
            tG.append(G[i] % p0)
        t1 = poly_mul(logn, tf, tG, p0)
        t2 = poly_mul(logn, tg, tF, p0)
        t3 = poly_sub(logn, t1, t2, p0)
        if t3[0] != q:
            return False
        for i in range(1, n):
            if t3[i] != 0:
                return False

    return (F, G)

# ==========================================================================
# Support function: SHAKE-based PRNG.
#
# This PRNG implementation uses SHAKE256 to generate some arbitrary-length
# data. The Python3 standard hashlib implementation of SHAKE256 does not
# support a streamable output; this code therefore regenerates its output
# repeatedly, and becomes quite inefficient over time. This is good enough
# for testing purposes, though; if generating many keys successively, then
# a new PRNG instance should be created for each keygen.

class SHAKE256x4:
    def __init__(self, seed):
        # We initialize 4 parallel SHAKE256 instances.
        # Each instance is fed with 'seed || i' with 'i' being the
        # instance number (0 to 3, over a single byte).
        sh = []
        for i in range(0, 4):
            s = hashlib.shake_256()
            s.update(seed)
            s.update(int(i).to_bytes(1, byteorder='little'))
            sh.append(s)
        self.sh = sh
        self.buf = []
        self.ptr = 0

    # Get the next 512-byte output block.
    def next_512(self):
        # We generate the stream up to some length, and buffer the output.
        # Owing to the limitations of hashlib's API, we need to regenerate
        # the whole stream whenever we have to obtain data beyond the
        # buffered output.
        if self.ptr == len(self.buf):
            kx = len(self.buf) << 1
            if kx == 0:
                kx = 512
            k1 = kx >> 2
            bb = []
            for i in range(0, 4):
                bb.append(self.sh[i].digest(k1))
            # The four SHAKE256 outputs are interleaved with 8-byte
            # granularity. This matches the C code; this interleaving
            # mimics what happens with running the four SHAKE256
            # simultaneously with AVX2 opcodes.
            self.buf = b''
            for i in range(0, kx >> 5):
                for j in range(0, 4):
                    self.buf += bb[j][(i << 3):((i << 3) + 8)]
        i = self.ptr
        j = i + 512
        self.ptr = j
        return self.buf[i:j]

# ==========================================================================
# Support functions: polynomials.
#
# All polynomials are considered modulo x^n+1, for some integer n which
# is a power of two  between 1 and 1024 (inclusive). Coefficients are
# integers, but they are implicitly reduced modulo a given big integer
# m, with signed convention: the output is in [-m/2..+m/2[ (if m is
# even, the value can be -m/2 but not +m/2). The modulus m depends on
# the usage context; it is either a power of two (2^(31*k) for some
# integer k), or a product of small primes p_i (primes lower than 2^31,
# with p_i = 1 mod 2048).
#
# The functions below simply use big integers and explicit modular
# reductions, for clarity of exposition. When m = 2^(31*k), each integer
# is assumed to consist of k limbs, each being an unsigned 31-bit integer
# (0 to 2^31-1); the least significant limb is limb 0, and the most
# significant is limb k-1. The top bit of limb k-1 (i.e. bit 30) is the
# "sign bit", using two's complement interpretation for negative values.
#
# When m is the product of k primes p_i, then a practical, constant-time
# implementation may switch to RNS, i.e. represent an integer z with
# the k-length sequence of values (z mod p_i), and perform computations
# modulo each p_i in parallel. This naturally implements the semantics
# of operations modulo m, and promotes efficiency since the p_i are
# small enough to fit in one CPU register (with an extra free bit, which
# is convenient for some operations). Moreover, since p_i = 1 mod 2048,
# the NTT can be used to further speed up some polynomial operations
# modulo p_i.
#
# The exact implementation method is not enforced, as long as it returns the
# correct results, i.e. with the reduction modulo m.
#
# Unless stated explicitly otherwise, all functions below assume that their
# inputs are already in the proper range for the intended modulus m.

# Reduce an integer modulo m, with signed convention.
def mred(z, m):
    z %= m
    hm = (m + 1) >> 1
    if z >= hm:
        z -= m
    return z

# Split a polynomial into its even and odd coefficients.
# If f = fe(x^2) + x*fo(x^2), then this function returns fe and fo, which
# are half-length polynomials.
def poly_split(logn, f, m):
    n = 1 << logn
    hn = n >> 1
    fe = []
    fo = []
    for i in range(0, hn):
        fe.append(f[2 * i])
        fo.append(f[2 * i + 1])
    return (fe, fo)

# Given polynomial f modulo x^n + 1, return f(x^2) modulo x^(2n) + 1.
def poly_expand(logn, f, m):
    n = 1 << logn
    h = [];
    for i in range(0, n):
        h.append(f[i])
        h.append(0)
    return h

# Add two polynomials modulo x^n + 1, with resulting coefficients
# normalized modulo the given integer m. The normalization is signed:
# coefficients are returned as integers in the [-m/2..+m/2[ range.
def poly_add(logn, f, g, m):
    n = 1 << logn
    h = []
    for i in range(0, n):
        h.append(mred(f[i] + g[i], m))
    return h

# Subtract two polynomials modulo x^n + 1, with resulting coefficients
# normalized modulo the given integer m. The normalization is signed:
# coefficients are returned as integers in the [-m/2..+m/2[ range.
def poly_sub(logn, f, g, m):
    n = 1 << logn
    h = []
    for i in range(0, n):
        h.append(mred(f[i] - g[i], m))
    return h

# Multiply two polynomials modulo x^n + 1, with resulting coefficients
# normalized modulo the given integer m. The normalization is signed:
# coefficients are returned as integers in the [-m/2..+m/2[ range.
def poly_mul(logn, f, g, m):
    # This uses a basic quadratic algorithm, which is of course slow.
    # Karatsuba-Ofman could be used, or, in practice, RNS+NTT (when
    # m is a product of small primes p_i).
    n = 1 << logn
    fg = []
    for k in range(0, 2*n):
        fg.append(0)
    for i in range(0, n):
        for j in range(0, n):
            fg[i + j] += f[i] * g[j]
    h = []
    for i in range(0, n):
        h.append(mred(fg[i] - fg[i + n], m))
    return h

# Multiply a polynomial by x (modulo x^n + 1).
def poly_mulx(logn, f, m):
    n = 1 << logn
    h = []
    hm = m >> 1
    if (m & 1) == 0 and f[n - 1] == -hm:
        # If m is even and the coefficient is -m/2, then negation
        # overflows and the coefficient is unchanged.
        h.append(f[n - 1])
    else:
        h.append(-f[n - 1])
    for i in range(0, n - 1):
        h.append(f[i])
    return h

# Get the Galois conjugate of a polynomial f, modulo x^n + 1. This is
# simply negation of odd-indexed coefficients.
def poly_conj(logn, f, m):
    n = 1 << logn
    h = []
    if (m & 1) == 0:
        # If m is even and the coefficient is -m/2, then negation
        # overflows and the coefficient is unchanged.
        hm = m >> 1
        for i in range(0, n, 2):
            h.append(f[i])
            if f[i + 1] == -hm:
                h.append(f[i + 1])
            else:
                h.append(-f[i + 1])
    else:
        for i in range(0, n, 2):
            h.append(f[i])
            h.append(-f[i + 1])
    return h

# Get the field norm of a polynomial f. Given f modulo x^n + 1, we
# compute N(f) such that N(f)(x^2) = f*conj(f). The Galois conjugate is
# obtained by negating the odd-indexed coefficients. In other words,
# given f, we first split it into its even-indexed and odd-indexed
# coefficients, as two half-size polynomials:
#    f = fe(x^2) + x*fo(x^2)
# Then:
#    conj(f) = fe(x^2) - x*fo(x^2)
#    N(f) = fe^2 - x*fo^2
# The coefficients of N(f) are normalized (signed) modulo the provided m.
def poly_field_norm(logn, f, m):
    (fe, fo) = poly_split(logn, f, m)
    he = poly_mul(logn - 1, fe, fe, m)
    ho = poly_mul(logn - 1, fo, fo, m)
    xho = poly_mulx(logn - 1, ho, m)
    return poly_sub(logn - 1, he, xho, m)

# Get the Hermitian adjoint of a polynomial f. If:
#   f = \sum_{i=0}^{n-1} f[i]*x^i
# Then its adjoint is:
#   adj(f) = f[0] - \sum_{i=1}^{n-1} f[n-i]*x^i
# In FFT representation this could be computed by complex conjugation of
# the complex coefficients. This function works with normal representation
# and integers.
def poly_adj(logn, f, m):
    n = 1 << logn
    h = []
    h.append(f[0])
    if (m & 1) == 0:
        # If m is even and the coefficient is -m/2, then negation
        # overflows and the coefficient is unchanged.
        hm = m >> 1
        for i in range(1, n):
            if f[n - i] == -hm:
                h.append(f[n - i])
            else:
                h.append(-f[n - i])
    else:
        for i in range(1, n):
            h.append(-f[n - i])
    return h

# Get some limbs out of an integer z. For the input integer z, this
# function returns floor((z mod 2^(31*(blen+rlen))) / 2^(31*blen))
# (this value is always between 0 and 2^(31*rlen) - 1).
# If rlen == 0 then this function returns 0.
def extract_limbs(z, blen, rlen):
    return (z >> (31*blen)) & ((1 << (31*rlen)) - 1)

# Test whether the top bit of a given unsigned integer (of size rlen
# limbs) is 1.
def is_negative(z, rlen):
    if rlen == 0:
        return False
    else:
        return (z >> (31*rlen - 1)) != 0

# Extract the maximum bit length of any coefficient in the provided
# polynomial. Only some 31-bit limbs are considered; namely:
#
#  - Each coefficient f[i] is considered to be represented as a
#    sequence of 31-bit limbs, with limb j being the value:
#       floor((f[i] mod 2^(31*(j + 1))) / (2^(31*j)))
#
#  - Only limbs j such that:
#       blen <= j < blen + rlen
#    are considered.
#
#  - If rlen == 0 then the function returns 0. Otherwise, if limb of
#    index blen+rlen-1 has its top bit set (i.e. the limb value is not lower
#    than 2^30), then all the limb bits are complemented (the provided
#    polynomial f is not modified; the complementation is local to this
#    function).
#
#  - After the potential complementation, the retained limbs encode an
#    unsigned integer lower than 2^(31*rlen - 1). The bit length of the
#    coefficient is the minimal size, in bits, of the binary representation
#    of this integer; it is thus between 0 and 31*rlen - 1 (inclusive).
#
# Use of complementation corresponds to a notion of bit length which is
# the minimal length of the binary representation of the integer, using
# two's complement for negative value, and _excluding_ the sign bit. In
# particular, for an integer k, 2^k has length k+1 bits, but -2^k has
# length k bits.
def poly_max_bitlength(logn, f, blen, rlen):
    if rlen == 0:
        return 0
    t = 0
    for i in range(0, 1 << logn):
        r = extract_limbs(f[i], blen, rlen)
        if is_negative(r, rlen):
            # For a k-bit integer z, inverting all bits is equivalent
            # to subtracting z from (2^k)-1.
            r = ((1 << (31*rlen)) - 1) - r
        t = max(t, int(r).bit_length())
    return t

# Convert an integer polynomial to fixed-point, with scaling. For
# each coefficient f[i]:
#   - We only consider limbs blen to blen+rlen-1 as part of the integer;
#     the top bit of limb blen+rlen-1 is used as the sign.
#   - If the extracted limbs, after sign processing, encode integer y,
#     then the corresponding coefficient in the output is an approximation
#     of y/2^sc.
def poly_to_fxr(logn, f, blen, rlen, sc):
    # For a practical implementation: blen and rlen are public, but sc is
    # secret. To achieve constant-time processing, the code should scan all
    # limbs of each coefficient from indices blen to blen+rlen-1.
    h = []
    if rlen == 0:
        for i in range(0, 1 << logn):
            h.append(0)
    else:
        for i in range(0, 1 << logn):
            r = extract_limbs(f[i], blen, rlen)
            if is_negative(r, rlen):
                r -= 1 << (31*rlen)
            y = (r << 32) >> sc
            h.append(fxr_from_scaled32(y))
    return h

# ==========================================================================
# Support function: fixed-point integers and FFT.
#
# Each fixed point value is represented by a 64-bit unsigned integer
# (value in [0..2^64-1]). The unsigned integer value u represents the
# real number v/2^32, where v is the signed 64-bit interpretation of u
# (using two's complement for negative values).

# Convert an integer to fixed-point. This tolerates (and truncates)
# out-of-range integers.
def fxr_from_int(t):
    return (t << 32) & 0xFFFFFFFFFFFFFFFF

# Round a fixed-point value to the nearest integer.
def fxr_round32(z):
    z = (z + 0x80000000) & 0xFFFFFFFFFFFFFFFF
    return fxr_to_scaled32(z) >> 32

# For integer z, return the fixed-point representation of z/2^32 (possibly
# truncated).
def fxr_from_scaled32(z):
    return z & 0xFFFFFFFFFFFFFFFF

# Convert a fixed-point value z to integer v such that v = z/2^32. This
# really is a sign extension.
def fxr_to_scaled32(z):
    return z - ((z & 0x8000000000000000) << 1)

# Addition of fixed-point values.
def fxr_add(z1, z2):
    return (z1 + z2) & 0xFFFFFFFFFFFFFFFF

# Subtraction of fixed-point values.
def fxr_sub(z1, z2):
    return (z1 - z2) & 0xFFFFFFFFFFFFFFFF

# Negation of a fixed-point value.
def fxr_neg(z):
    return (-z) & 0xFFFFFFFFFFFFFFFF

# Get the absolute value of a fixed-point value.
def fxr_abs(z):
    nm = (z >> 63) & 1
    return (z - ((z << 1) & -nm)) & 0xFFFFFFFFFFFFFFFF

# Multiplication of fixed-point values.
def fxr_mul(z1, z2):
    i1 = fxr_to_scaled32(z1)
    i2 = fxr_to_scaled32(z2)
    return ((i1 * i2) >> 32) & 0xFFFFFFFFFFFFFFFF

# Squaring of a fixed-point value.
def fxr_sqr(z):
    return fxr_mul(z, z)

# Division of a fixed-point value by 2. Rounding is applied.
def fxr_half(z):
    z = (z + 1) & 0xFFFFFFFFFFFFFFFF
    i = fxr_to_scaled32(z)
    return (i >> 1) & 0xFFFFFFFFFFFFFFFF

# Multiplication of a fixed-point value by a power of 2. The value is
# multiplied by 2^e. Exponent e shall be in the 0 to 30 range.
def fxr_mul2e(z, e):
    return (z << e) & 0xFFFFFFFFFFFFFFFF

# Inversion of a fixed-point value.
def fxr_inv(z):
    return fxr_div(1 << 32, z)

# Fixed-point division.
def fxr_div(z1, z2):
    # We mimic the C code, so that we get the exact same output in all cases.
    # In any case, a practical implementation would need to be constant-time
    # and thus avoid the division opcode of the CPU.
    s1 = z1 >> 63
    s2 = z2 >> 63
    z1 = fxr_abs(z1)
    z2 = fxr_abs(z2)
    qq = 0
    num = z1 >> 31
    for i in range(63, 32, -1):
        b = 1 - (((num - z2) & 0xFFFFFFFFFFFFFFFF) >> 63)
        qq |= b << i
        num -= z2 & -b
        num <<= 1
        num |= (z1 >> (i - 33)) & 1
    for i in range(32, -1, -1):
        b = 1 - (((num - z2) & 0xFFFFFFFFFFFFFFFF) >> 63)
        qq |= b << i
        num -= z2 & -b
        num <<= 1

    # Rounding
    b = 1 - (((num - z2) & 0xFFFFFFFFFFFFFFFF) >> 63)
    qq += b

    # Sign management
    qq -= (qq << 1) & -((s1 + s2) & 1)
    return qq & 0xFFFFFFFFFFFFFFFF

# We represent complex numbers as pairs of fixed-point values (real and
# imaginary parts).

# Fixed-point complex addition.
def fxc_add(c1, c2):
    (r1, i1) = c1
    (r2, i2) = c2
    return (fxr_add(r1, r2), fxr_add(i1, i2))

# Fixed-point complex subtraction.
def fxc_sub(c1, c2):
    (r1, i1) = c1
    (r2, i2) = c2
    return (fxr_sub(r1, r2), fxr_sub(i1, i2))

# Fixed-point complex halving.
def fxc_half(c):
    (r, i) = c
    return (fxr_half(r), fxr_half(i))

# Fixed-point complex multiplication.
def fxc_mul(c1, c2):
    (r1, i1) = c1
    (r2, i2) = c2
    # As in the C code, we use the 3-mul formulas.
    t0 = fxr_mul(r1, r2)
    t1 = fxr_mul(i1, i2)
    t2 = fxr_mul(fxr_add(r1, i1), fxr_add(r2, i2))
    return (fxr_sub(t0, t1), fxr_sub(t2, fxr_add(t0, t1)))

# Fixed-point complex conjugation (i.e. negation of the imaginary part).
def fxc_conj(c):
    (r, i) = c
    return (r, fxr_neg(i))

# We need a table of 2048-th roots of unity, in the proper order for the
# FFT, i.e. with "bit reversal". For i = 0 to 1023, gm_tab[i] contains
# exp(rev10(i)*pi/1024), as a complex number in our fixed-point
# representation. rev10() is the bit-reversal function over 10 bits.
gm_tab = []
def mk_FFT_root():
    # We hash generated values for consistency check.
    # The code below uses math.cos() and math.sin(), which work with
    # the usual IEEE754 'binary64' type (precision is 53 bits); this
    # is enough to obtain the correct values (the reference hash was
    # obtained from Sage computations with enough precision to avoid
    # any ambiguous rounding).
    h = hashlib.sha256()
    for i in range(0, 1024):
        j = 0
        for k in range(0, 10):
            j += ((i >> k) & 1) << (9 - k)
        x = round(math.cos(j*math.pi/1024)*2**32) & 0xFFFFFFFFFFFFFFFF
        y = round(math.sin(j*math.pi/1024)*2**32) & 0xFFFFFFFFFFFFFFFF
        gm_tab.append((x, y))
        h.update(int(x).to_bytes(8, byteorder='little'))
        h.update(int(y).to_bytes(8, byteorder='little'))
    assert h.hexdigest() == 'e0746da83816f05bd5a3bbe0867bba11ee2bab9cd781780ce4d205d6c51b03d9'

mk_FFT_root()

# Convert a polynomial into FFT representation. The source polynomial
# has n coefficients, which are all real numbers (in fixed-point
# representation). In FFT representation, it has n/2 complex
# coefficients, and coefficient j has its real value in index j and its
# imaginary value in index j+n/2. Only n/2 coefficients are used because
# for a source real polynomial, because in that case the full FFT (n
# complex coefficients) is redundant: n/2 of the coefficients are
# conjugates of the other n/2. The separation of the real and imaginary
# part allows optimizing out the first iteration of the FFT loop.
#
# The provided polynomial is converted in place.
def vect_FFT(logn, f):
    hn = 1 << (logn - 1)
    t = hn
    for lm in range(1, logn):
        m = 1 << lm
        ht = t >> 1
        j0 = 0
        hm = m >> 1
        for i in range(0, hm):
            s = gm_tab[m + i]
            for j in range(j0, j0 + ht):
                x = (f[j], f[j + hn])
                y = (f[j + ht], f[j + ht + hn])
                z = fxc_mul(s, y)
                (f[j], f[j + hn]) = fxc_add(x, z)
                (f[j + ht], f[j + ht + hn]) = fxc_sub(x, z)
            j0 += t
        t = ht

# Inverse FFT: does the opposite operation of fx_FFT()
# (since fixed-point values imply rounding, using fx_FFT() then fx_iFFT()
# does not, in general, yield back the _exact_ same polynomial as the
# original input to fx_FFT(), only a close approximation).
#
# The provided polynomial is converted in place.
#
# Note: in the final outer iteration:
#
#  - f[j] and f[j+hn] are set to the half of a complex number. Fixed point
#    values are held in integers in the [-2^63..+2^63-1] range (in signed
#    interpretation) and halving is done by adding 1 and then an
#    arithmetic right shift, thus with an output necessarily in the
#    [-2^62..+2^62-1] range (addition of 1 to 2^63-1 "wraps around" in that
#    +2^63 is interpreted as -2^63 instead).
#
#  - f[j+ht] and f[j+ht+hn] are set to the product of the half of a
#    complex number, and the complex sqrt(2)-i*sqrt(2). Value sqrt(2)
#    is represented by the fixed point 3037000500. Following the steps
#    in fxc_mul():
#       r1, i1:           -2^62 .. +2^62-1
#       r2:               +3037000500
#       i2:               -3037000500
#       t0:               -3260954456358912000 .. +3260954456358911999
#       t1:               -3260954456358912000 .. +3260954456358912000
#       fxr_add(r1, i1):  -2^63 .. +2^63-2
#       fxr_add(r2, i2):  always zero
#       t2:               always zero
#       fxr_sub(t0, t1):  -6521908912717824000 .. +6521908912717823999
#       fxr_add(t0, t1):  -6521908912717824000 .. +6521908912717823999
#    Thus, the obtained output values must be in the
#    [-6521908912717824000..+6521908912717824000] range.
#
# If the output of vect_iFFT() is then rounded to integers, then the
# maximum range for any output value after rounding is
# [-1518500250..+1518500250], regardless of the contents of the source
# vector. In particular, this is a smaller range than signed 32-bit integers
# in general, and it avoids the troublesome values such as -2^31.
def vect_iFFT(logn, f):
    hn = 1 << (logn - 1)
    ht = 1
    for lm in range(logn - 1, 0, -1):
        m = 1 << lm
        t = ht << 1
        j0 = 0
        hm = m >> 1
        for i in range(0, hm):
            s = fxc_conj(gm_tab[m + i])
            for j in range(j0, j0 + ht):
                x = (f[j], f[j + hn])
                y = (f[j + ht], f[j + ht + hn])
                z = fxc_half(fxc_sub(x, y))
                (f[j], f[j + hn]) = fxc_half(fxc_add(x, y))
                (f[j + ht], f[j + ht + hn]) = fxc_mul(s, z)
            j0 += t
        ht = t

# Convert an integer polynomial (with small coefficients) into a real
# polynomial (with fixed-point representation).
def vect_to_fxr(logn, f):
    v = []
    for i in range(0, 1 << logn):
        v.append(fxr_from_int(f[i]))
    return v

# Add two real polynomials together.
def vect_add(logn, f, g):
    h = []
    for i in range(0, 1 << logn):
        h.append(fxr_add(f[i], g[i]))
    return h

# Multiply a polynomial with a real constant (c).
def vect_mul_realconst(logn, f, c):
    h = []
    for i in range(0, 1 << logn):
        h.append(fxr_mul(f[i], c))
    return h

# Multiply a polynomial with a power of two. Exponent e should be in
# the 0 to 30 range.
def vect_mul2e(logn, f, e):
    h = []
    for i in range(0, 1 << logn):
        h.append(fxr_mul2e(f[i], e))
    return h

# Multiply two polynomials together. The source polynomials must be in FFT
# representation, and the output is in FFT representation.
def vect_mul_fft(logn, f, g):
    hn = 1 << (logn - 1)
    h = [0]*(hn << 1)
    for i in range(0, hn):
        (h[i], h[i + hn]) = fxc_mul((f[i], f[i + hn]), (g[i], g[i + hn]))
    return h

# Convert a polynomial (in FFT representation) to its Hermitian adjoint
# (i.e. each coefficient is replaced with its conjugate).
#
# The provided polynomial is converted in place.
def vect_adj_fft(logn, f):
    hn = 1 << (logn - 1)
    h = []
    for i in range(0, hn):
        h.append(f[i])
    for i in range(hn, hn << 1):
        h.append(fxr_neg(f[i]))
    return h

# Multiply two polynomials f and g, both in FFT representation; the
# second polynomial (g) is assumed to be self-adjoint, i.e. equal to
# its adjoint (this implies that all its FFT coefficients are real).
# Array g[] may have half-length.
def vect_mul_selfadj_fft(logn, f, g):
    hn = 1 << (logn - 1)
    h = []
    for i in range(0, hn):
        h.append(fxr_mul(f[i], g[i]))
    for i in range(hn, hn << 1):
        h.append(fxr_mul(f[i], g[i - hn]))
    return h

# Divide polynomial f by polynomial g, both in FFT representation; the
# second polynomial (g) is assumed to be self-adjoint, i.e. equal to
# its adjoint (this implies that all its FFT coefficients are real).
# Array g[] may have half-length.
def vect_div_selfadj_fft(logn, f, g):
    hn = 1 << (logn - 1)
    h = []
    for i in range(0, hn):
        h.append(fxr_div(f[i], g[i]))
    for i in range(hn, hn << 1):
        h.append(fxr_div(f[i], g[i - hn]))
    return h

# Compute f*adj(f) + g*adj(g). Input and output polynomials are in FFT
# representation. Output polynomial is self-adjoint and is returned as
# an half-length polynomial (the imaginary parts of the output coefficients
# are implicitly all equal to zero).
def vect_norm_fft(logn, f, g):
    hn = 1 << (logn - 1)
    h = []
    for i in range(0, hn):
        z1 = fxr_add(fxr_sqr(f[i]), fxr_sqr(f[i + hn]))
        z2 = fxr_add(fxr_sqr(g[i]), fxr_sqr(g[i + hn]))
        h.append(fxr_add(z1, z2))
    return h

# Compute 2^e/(f*adj(f) + g*adj(g)). Input and output polynomials are in FFT
# representation. Output polynomial is self-adjoint and is returned as
# an half-length polynomial (the imaginary parts of the output coefficients
# are all equal to zero).
def vect_invnorm_fft(logn, f, g, e):
    hn = 1 << (logn - 1)
    h = []
    r = fxr_from_int(1 << e)
    for i in range(0, hn):
        z1 = fxr_add(fxr_sqr(f[i]), fxr_sqr(f[i + hn]))
        z2 = fxr_add(fxr_sqr(g[i]), fxr_sqr(g[i + hn]))
        h.append(fxr_div(r, fxr_add(z1, z2)))
    return h

# ==========================================================================
# Support function: NTT modulo x^n+1 and modulo q = 12289.
#
# This code is used to verify that the public key h = g/f mod x^n+1 mod q
# is computable, and also to compute it.

# For g = 7 (which is a primitive 2048-th root of 1 modulo q = 12289),
# we make two precomputed tables:
#   mq_gm[i] = g^(rev10(i)) mod q
#   mq_igm[i] = 1/(2*g^(rev10(i))) mod q
mq_gm = [0] * 1024
mq_igm = [0] * 1024

def mk_mqgm():
    n = 1 << 10
    g = 7      # primitive 2048-th root of 1 modulo q
    ig = 8778  # 1/g mod q
    x1 = 1
    x2 = (q + 1) >> 1
    for i in range(0, n):
        # j = rev10(i)
        j = (i >> 5) | (i << 5)
        jm = j & 0b0010000100
        j = ((j & 0b1100011000) >> 3) | ((j & 0b0001100011) << 3)
        j = ((j & 0b1001010010) >> 1) | ((j & 0b0100101001) << 1) | jm
        mq_gm[j] = x1
        mq_igm[j] = x2
        x1 = (g * x1) % q
        x2 = (ig * x2) % q

mk_mqgm()

# Apply NTT modulo x^n + 1 and modulo q = 12289.
def mq_NTT(logn, f):
    assert 0 <= logn and logn <= 10
    assert len(f) == 1 << logn

    # Nothing to do if degree is 1.
    if logn == 0:
        return

    # logn outer iterations.
    t = 1 << logn
    for lm in range(0, logn):
        m = 1 << lm
        ht = t >> 1
        v0 = 0
        for u in range(0, m):
            s = mq_gm[u + m]
            for v in range(0, ht):
                k1 = v0 + v
                k2 = k1 + ht
                x1 = f[k1]
                x2 = f[k2] * s
                f[k1] = (x1 + x2) % q
                f[k2] = (x1 - x2) % q
            v0 += t
        t = ht

# Opposite transform of mq_NTT().
def mq_iNTT(logn, f):
    assert 0 <= logn and logn <= 10
    assert len(f) == 1 << logn

    # Nothing to do if degree is 1.
    if logn == 0:
        return

    # logn outer iterations, just reversing the work of mq_NTT().
    # Note: mq_igm[] already contains a 1/2 factor baked in.
    t = 1
    for lm in range(0, logn):
        hm = 1 << (logn - 1 - lm)
        dt = t << 1
        v0 = 0
        for u in range(0, hm):
            s = mq_igm[u + hm]
            for v in range(0, t):
                k1 = v0 + v
                k2 = k1 + t
                x1 = f[k1]
                x2 = f[k2]
                z = (x1 + x2) % q
                z += q & -(z & 1)
                f[k1] = z >> 1
                f[k2] = ((x1 - x2) * s) % q
            v0 += dt
        t = dt

# Division modulo q = 12289.
# We do not have to support the case of a zero divisor (caller ensures that
# this case does not happen).
def mq_div(x, y):
    # Using Fermat's Little Theorem: 1/y = y^(q-2) mod q.
    y2 = (y * y) % q
    y3 = (y * y2) % q
    y5 = (y2 * y3) % q
    y10 = (y5 * y5) % q
    y20 = (y10 * y10) % q
    y40 = (y20 * y20) % q
    y80 = (y40 * y40) % q
    y160 = (y80 * y80) % q
    y163 = (y160 * y3) % q
    y323 = (y163 * y160) % q
    y646 = (y323 * y323) % q
    y1292 = (y646 * y646) % q
    y1455 = (y1292 * y163) % q
    y2910 = (y1455 * y1455) % q
    y5820 = (y2910 * y2910) % q
    y6143 = (y5820 * y323) % q
    y12286 = (y6143 * y6143) % q
    y12287 = (y12286 * y) % q
    return (x * y12287) % q

# ==========================================================================

# A test function. It generates 100 key pairs for the specified degree,
# with a given seed ('test0' to 'test99'). For each key pair, it
# verifies that h is the correct public key for f and g. Moreover, the
# obtained private key (f, g, F, G) is hashed with SHAKE256 and the
# result is output in hexadecimal, for comparing with keys obtained with
# the C implementation.
def _run_test(logn):
    for i in range(0, 100):
        n = 1 << logn
        seed = bytes('test%d' % i, 'utf-8')
        prng = SHAKE256x4(seed)
        (f, g, F, G, h) = Falcon_keygen(logn, prng)

        # Check that f*h - g = 0 mod x^n+1 mod q.
        for j in range(0, n):
            assert 0 <= h[j] and h[j] < q
        fh = [0] * (2*n)
        for j in range(0, n):
            for k in range(0, n):
                fh[j + k] += f[j] * h[k]
        for j in range(0, n):
            z = (fh[j] - fh[j + n] - g[j]) % q
            assert z == 0

        # Hash and print the private key.
        bf = b''
        bg = b''
        bF = b''
        bG = b''
        for j in range(0, n):
            bf += int(f[j] & 0xFF).to_bytes(1, byteorder='little')
            bg += int(g[j] & 0xFF).to_bytes(1, byteorder='little')
            bF += int(F[j] & 0xFF).to_bytes(1, byteorder='little')
            bG += int(G[j] & 0xFF).to_bytes(1, byteorder='little')
        sh = hashlib.shake_256()
        sh.update(bf)
        sh.update(bg)
        sh.update(bF)
        sh.update(bG)
        print('%4d: %s' % (i, sh.hexdigest(20)), flush=True)

import sys
if __name__ == '__main__':
    aa = sys.argv
    if len(aa) < 2:
        logn = 9
    else:
        logn = int(aa[1])
    _run_test(logn)
