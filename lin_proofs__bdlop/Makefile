CPP = g++ -g3
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -I NFLlib/include/ -I NFLlib/include/nfl -I NFLlib/include/nfl/prng -I include -DNFL_OPTIMIZED=ON -DNTT_AVX2 -flto=auto
INCLUDES = include/bench.h include/cpucycles.h
BLAKE3 = src/blake3/blake3.c src/blake3/blake3_dispatch.c src/blake3/blake3_portable.c src/blake3/blake3_sse2_x86-64_unix.S src/blake3/blake3_sse41_x86-64_unix.S src/blake3/blake3_avx2_x86-64_unix.S src/blake3/blake3_avx512_x86-64_unix.S
BENCH = src/bench.c src/cpucycles.c
TEST = src/test.c
LIBS = deps/libnfllib_static.a -lgmp -lmpfr -L deps/ -lflint -lquadmath

all: olingo

olingo: src/olingo.cpp src/bdlop.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/bdlop.cpp -o bdlop_olingo.o
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small_olingo.o
	${CPP} ${CFLAGS} -DMAIN src/olingo.cpp bdlop_olingo.o sample_z_small_olingo.o ${TEST} ${BENCH} ${BLAKE3} -o olingo ${LIBS}

clean:
	rm *.o olingo