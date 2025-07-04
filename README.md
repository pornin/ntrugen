# ntrugen

This is an implementation of the key pair generation process used in
three lattice-based cryptographic schemes that use NTRU lattices, namely
[BAT](https://eprint.iacr.org/2022/031),
[Falcon](https://falcon-sign.info/), and
[Hawk](https://hawk-sign.info/). This implementation is meant
for research and integration into implementations of these schemes; it
has the following characteristics:

  - It is written in standard C (C99) and is portable. It can be run on
    large systems but also on small microcontrollers.
  - It includes numerous optional optimizations for specific platforms,
    in particular x86 CPUs with AVX2 intrinsics.
  - It is fully constant-time.
  - It uses no floating-point operation whatsoever.
  - It is faster and uses less RAM than previously published
    implementations.

See [the accompanying note](ntrugen.pdf) (also on
[eprint](https://eprint.iacr.org/2023/290)) and the [updated
note](ntrugen2.pdf) for details on the involved techniques, and some
benchmarks.

**Warning:** all the implemented algorithms may evolve in the near future;
in particular, Falcon is currently undergoing specification by NIST. I
will try to track these variants in the implementation, which may imply
some changes in the behaviour and API.
