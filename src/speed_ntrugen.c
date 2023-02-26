#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ntrugen.h"
#include "ntrugen_prng.h"

#ifndef DO_BENCH86
#if defined __i386__ || defined _M_IX86 || defined __x86_64__ || defined _M_X64
#define DO_BENCH86   1
#else
#define DO_BENCH86   0
#endif
#endif

#if DO_BENCH86
#include <immintrin.h>

#if defined __GNUC__ || defined __clang__
__attribute__((target("sse2")))
#endif
static inline uint64_t
core_cycles(void)
{
#if defined __GNUC__ && !defined __clang__
	uint32_t hi, lo;

	_mm_lfence();
	__asm__ __volatile__ ("rdtsc" : "=d" (hi), "=a" (lo) : : );
	return ((uint64_t)hi << 32) | (uint64_t)lo;
#else
	_mm_lfence();
	return __rdtsc();
#endif
}

#endif

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

/*
 * Benchmark function takes an opaque context and an iteration count;
 * it returns 0 on success, a negative error code on error.
 */
typedef int (*bench_fun)(void *ctx, unsigned long num);

/*
 * Returned value is the time per iteration in nanoseconds.
 * WARNING: ON x86, VALUES ARE RETURNED IN CLOCK CYCLES, NOT NANOSECONDS;
 * THRESHOLD IS IN BILLIONS OF CYCLES.
 *
 * If the benchmark function reports an error, 0.0 is returned.
 */
static double
do_bench(bench_fun bf, void *ctx, double threshold)
{
	unsigned long num;
	int r;

	/*
	 * Always do a few blank runs to "train" the caches and branch
	 * prediction.
	 */
	r = bf(ctx, 5);
	if (r != 0) {
		fprintf(stderr, "ERR: %d\n", r);
		return 0.0;
	}

	num = 1;
	for (;;) {
#if DO_BENCH86
		uint64_t begin, end;
#else
		clock_t begin, end;
#endif
		double tt;

#if DO_BENCH86
		begin = core_cycles();
#else
		begin = clock();
#endif
		r = bf(ctx, num);
#if DO_BENCH86
		end = core_cycles();
#else
		end = clock();
#endif
		if (r != 0) {
			fprintf(stderr, "ERR: %d\n", r);
			return 0.0;
		}
#if DO_BENCH86
		tt = (double)(end - begin) / (double)1000000000.0;
#else
		tt = (double)(end - begin) / (double)CLOCKS_PER_SEC;
#endif
		if (tt >= threshold) {
			return tt * 1000000000.0 / (double)num;
		}

		/*
		 * If the function ran for less than 0.3 seconds then
		 * we simply double the iteration number; otherwise, we
		 * use the run time to try to get a "correct" number of
		 * iterations quickly.
		 */
		if (tt < 0.3) {
			num <<= 1;
		} else {
			unsigned long num2;

			num2 = (unsigned long)((double)num
				* (threshold * 1.1) / tt);
			if (num2 <= num) {
				num2 = num + 1;
			}
			num = num2;
		}
		num <<= 1;
	}
}

typedef struct {
	uint64_t tmp[3 * 1024 + 8];
	int8_t f[1024], g[1024];
	ntrugen_prng_shake_context sc;
	unsigned logn;
	unsigned long last_num;
} bench_ntrugen_state;

static int
bench_ntrugen_BAT(void *ctx, unsigned long num)
{
	bench_ntrugen_state *bs = ctx;
	bs->last_num = num;
	while (num -- > 0) {
		int err = ntrugen_BAT_keygen(bs->logn, bs->f, bs->g,
			NULL, NULL, NULL,
			&ntrugen_prng_shake_out, &bs->sc,
			bs->tmp, sizeof bs->tmp);
		if (err < 0) {
			return err;
		}
	}
	return 0;
}

static int
bench_ntrugen_Falcon(void *ctx, unsigned long num)
{
	bench_ntrugen_state *bs = ctx;
	bs->last_num = num;
	while (num -- > 0) {
		int err = ntrugen_Falcon_keygen(bs->logn, bs->f, bs->g,
			NULL, NULL,
			&ntrugen_prng_shake_out, &bs->sc,
			bs->tmp, sizeof bs->tmp);
		if (err < 0) {
			return err;
		}
	}
	return 0;
}

static int
bench_ntrugen_Hawk(void *ctx, unsigned long num)
{
	bench_ntrugen_state *bs = ctx;
	bs->last_num = num;
	while (num -- > 0) {
		int err = ntrugen_Hawk_keygen(bs->logn, bs->f, bs->g,
			NULL, NULL, NULL, NULL, NULL, NULL,
			&ntrugen_prng_shake_out, &bs->sc,
			bs->tmp, sizeof bs->tmp);
		if (err < 0) {
			return err;
		}
	}
	return 0;
}

int
main(int argc, char *argv[])
{
	double threshold, scale;

	if (argc < 2) {
		threshold = 5.0;
	} else {
		threshold = atof(argv[1]);
	}
	if (threshold <= 0.0 || threshold > 60.0) {
		fprintf(stderr,
"usage: speed [ threshold ]\n"
"'threshold' is the minimum time for a bench run, in seconds (must be\n"
"positive and less than 60).\n");
		exit(EXIT_FAILURE);
	}
#if DO_BENCH86
	printf("time threshold = %.4f Gcyc\n", threshold);
#else
	printf("time threshold = %.4f s\n", threshold);
#endif
#if DO_BENCH86
	printf("x86 PLATFORM, USING TSC; VALUES IN MILLION CLOCK CYCLES\n");
	printf("(These values are MEANINGLESS if you did not disable"
		" frequency scaling.)\n");
	scale = 1 / 1000000.0;
#else
	printf("keygen in microseconds\n");
	scale = 1 / 1000.0;
#endif
	bench_ntrugen_state *bs = xmalloc(sizeof *bs);
	if (argc >= 3) {
		ntrugen_prng_shake_init(&bs->sc, argv[2], strlen(argv[2]));
	} else {
		uint8_t seed[32];
		if (ntrugen_sysrng(seed, sizeof seed) != 0) {
			fprintf(stderr, "system RNG failed\n");
			exit(EXIT_FAILURE);
		}
		ntrugen_prng_shake_init(&bs->sc, seed, sizeof seed);
	}

	double m;
	bs->logn = 8;
	m = do_bench(&bench_ntrugen_BAT, bs, threshold);
	printf("BAT-128-256     %9.3f  (%lu)\n", scale * m, bs->last_num);
	bs->logn = 9;
	m = do_bench(&bench_ntrugen_BAT, bs, threshold);
	printf("BAT-257-512     %9.3f  (%lu)\n", scale * m, bs->last_num);
	bs->logn = 10;
	m = do_bench(&bench_ntrugen_BAT, bs, threshold);
	printf("BAT-769-1024    %9.3f  (%lu)\n", scale * m, bs->last_num);

	bs->logn = 8;
	m = do_bench(&bench_ntrugen_Falcon, bs, threshold);
	printf("Falcon-256      %9.3f  (%lu)\n", scale * m, bs->last_num);
	bs->logn = 9;
	m = do_bench(&bench_ntrugen_Falcon, bs, threshold);
	printf("Falcon-512      %9.3f  (%lu)\n", scale * m, bs->last_num);
	bs->logn = 10;
	m = do_bench(&bench_ntrugen_Falcon, bs, threshold);
	printf("Falcon-1024     %9.3f  (%lu)\n", scale * m, bs->last_num);

	bs->logn = 8;
	m = do_bench(&bench_ntrugen_Hawk, bs, threshold);
	printf("Hawk-256        %9.3f  (%lu)\n", scale * m, bs->last_num);
	bs->logn = 9;
	m = do_bench(&bench_ntrugen_Hawk, bs, threshold);
	printf("Hawk-512        %9.3f  (%lu)\n", scale * m, bs->last_num);
	bs->logn = 10;
	m = do_bench(&bench_ntrugen_Hawk, bs, threshold);
	printf("Hawk-1024       %9.3f  (%lu)\n", scale * m, bs->last_num);

	xfree(bs);
	return 0;
}
