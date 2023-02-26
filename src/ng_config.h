/* ==================================================================== */
/*
 * Manual configuration.
 *
 * Each compile-time option is a macro that can defined to 1 (to enable
 * the feature) or 0 (to disable the feature). If the macro is not defined
 * then a default setting is applied (possibly with autodetection).
 */

#ifndef NG_CONFIG_H__
#define NG_CONFIG_H__

/*
 * OS-provided random source. Three sources are defined:
 *  - getentropy() system/library call
 *  - /dev/urandom special file
 *  - RtlGenRandom() on Windows
 * Default behaviour is to use getentropy() on Linux with Glibc-2.25+,
 * FreeBSD 12+ and OpenBSD; /dev/urandom on any other Unix-like system
 * (including Linux, FreeBSD, NetBSD, OpenBSD, DragonFly, macOS, Android,
 * Solaris, AIX); RtlGenRandom() on Windows.
 *
#define NTRUGEN_RAND_GETENTROPY   1
#define NTRUGEN_RAND_URANDOM      1
#define NTRUGEN_RAND_WIN32        1
 */

/*
 * AVX2 support. If defined to 1, then this enables use of AVX2 intrinsics,
 * even if the compiler does not explicitly target AVX2-only architectures.
 * Use of AVX2 intrinsics makes the compiled output not compatible with CPUs
 * that do not support AVX2.
 * Defining NTRUGEN_AVX2 to 0 forcibly disables use of AVX2 intrinsics
 * (the compiler might still use AVX2 opcodes as part of its automatic
 * vectorization of loops).
 * If NTRUGEN_AVX2 is not defined here, then use of AVX2 intrinsics will be
 * automatically detected based on the compilation target architecture
 * (e.g. use of the '-mavx2' compiler flag). There is no runtime detection
 * of the abilities of the CPU that actually runs the code.
 *
#define NTRUGEN_AVX2   1
 */

/*
 * ARM assembly support (for the ARMv7-M architecture with the optional
 * "DSP" extension, i.e. the ARM Cortex-M4 CPU). This assumes a GCC or
 * Clang compatible compiler.
 * If NTRUGEN_ASM_CORTEXM4 is not defined here, then it will be enabled
 * automatically if support is detected at compile-time (through
 * compiler-prodived predefined macros).
 *
#define NTRUGEN_ASM_CORTEXM4   1
 */

/*
 * Name prefixing. If this macro is defined to a identifier, then all
 * symbols corresponding to non-static functions and global data will
 * have a name starting with that identifier, followed by an underscore.
 * If undefined, then the default prefix is 'ntrugen'.
 * This is meant to be used for integration into larger applications, to
 * avoid name collisions. Also, the outer application might compile this
 * code twice, with and without AVX2 support, with different name prefixes,
 * and link both versions in the runtime, selecting at runtime which one to
 * call based on the current CPU abilities.
 *
#define NTRUGEN_PREFIX   ntrugen
 */

#endif
