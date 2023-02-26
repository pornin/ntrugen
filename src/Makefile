# Makefile for Unix-like systems (e.g. Linux).
#
# To use the AVX2 code, you can add the following to CFLAGS:
#   -DNTRUGEN_AVX2   enable the specific AVX2-aware code
#   -mavx2           blanket enabling of AVX2 for all code
#
# If -mavx2 is set, then this not only enables the AVX2-aware code
# (i.e. like setting -DNTRUGEN_AVX2), but also allows the compiler to
# use AVX2 opcodes for all other functions. With -DNTRUGEN_AVX2 only,
# only the internal functions that explicitly use AVX2 intrisics will
# use AVX2 opcodes, while the rest of the functions will keep complying
# with the standard ABI.

.POSIX:

CC = clang
CFLAGS = -Wall -Wextra -Wshadow -Wundef -O2 # -DNTRUGEN_AVX2
LD = clang
LDFLAGS =
LIBS =

OBJ = ng_bat.o ng_falcon.o ng_fxp.o ng_gauss.o ng_hawk.o ng_mp31.o ng_ntru.o ng_poly.o ng_zint31.o ng_prng.o sha3.o
TESTOBJ = test_ntrugen.o
SPEEDOBJ = speed_ntrugen.o

all: test_ntrugen speed_ntrugen

clean:
	-rm -f $(OBJ) $(TESTOBJ) test_ntrugen $(SPEEDOBJ) speed_ntrugen

test_ntrugen: $(OBJ) $(TESTOBJ)
	$(LD) $(LDFLAGS) -o test_ntrugen $(OBJ) $(TESTOBJ) $(LIBS)

speed_ntrugen: $(OBJ) $(SPEEDOBJ)
	$(LD) $(LDFLAGS) -o speed_ntrugen $(OBJ) $(SPEEDOBJ) $(LIBS)

ng_bat.o: ng_bat.c ng_inner.h ng_config.h sha3.h
	$(CC) $(CFLAGS) -c -o ng_bat.o ng_bat.c

ng_falcon.o: ng_falcon.c ng_inner.h ng_config.h sha3.h
	$(CC) $(CFLAGS) -c -o ng_falcon.o ng_falcon.c

ng_fxp.o: ng_fxp.c ng_inner.h ng_config.h sha3.h
	$(CC) $(CFLAGS) -c -o ng_fxp.o ng_fxp.c

ng_gauss.o: ng_gauss.c ng_inner.h ng_config.h sha3.h
	$(CC) $(CFLAGS) -c -o ng_gauss.o ng_gauss.c

ng_hawk.o: ng_hawk.c ng_inner.h ng_config.h sha3.h
	$(CC) $(CFLAGS) -c -o ng_hawk.o ng_hawk.c

ng_mp31.o: ng_mp31.c ng_inner.h ng_config.h sha3.h
	$(CC) $(CFLAGS) -c -o ng_mp31.o ng_mp31.c

ng_ntru.o: ng_ntru.c ng_inner.h ng_config.h sha3.h
	$(CC) $(CFLAGS) -c -o ng_ntru.o ng_ntru.c

ng_poly.o: ng_poly.c ng_inner.h ng_config.h sha3.h
	$(CC) $(CFLAGS) -c -o ng_poly.o ng_poly.c

ng_zint31.o: ng_zint31.c ng_inner.h ng_config.h sha3.h
	$(CC) $(CFLAGS) -c -o ng_zint31.o ng_zint31.c

sha3.o: sha3.c sha3.h
	$(CC) $(CFLAGS) -c -o sha3.o sha3.c

ng_prng.o: ng_prng.c ng_inner.h ng_config.h ntrugen.h ntrugen_prng.h sha3.h
	$(CC) $(CFLAGS) -c -o ng_prng.o ng_prng.c

test_ntrugen.o: test_ntrugen.c ng_inner.h ng_config.h ntrugen.h ntrugen_prng.h sha3.h
	$(CC) $(CFLAGS) -c -o test_ntrugen.o test_ntrugen.c

speed_ntrugen.o: speed_ntrugen.c ntrugen.h ntrugen_prng.h
	$(CC) $(CFLAGS) -c -o speed_ntrugen.o speed_ntrugen.c
