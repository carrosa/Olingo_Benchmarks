CC ?= /usr/bin/cc
THRESHOLD ?=
USERS ?=
DEFINES ?= -DTHRESHOLD=$(THRESHOLD) -DUSERS=$(USERS)
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wvla -Wpointer-arith -O3 -fomit-frame-pointer $(DEFINES)
NISTFLAGS += -Wno-unused-result -O3 -fomit-frame-pointer
SOURCES = sign.c packing.c polyvec.c poly.c ntt.c reduce.c rounding.c
HEADERS = config.h params.h api.h sign.h packing.h polyvec.h poly.h ntt.h \
  reduce.h rounding.h symmetric.h randombytes.h
KECCAK_SOURCES = $(SOURCES) fips202.c symmetric-shake.c
KECCAK_HEADERS = $(HEADERS) fips202.h
SOURCES256 = params.c poly256.c ntt256.c reduce256.c fastrandombytes.c test/cpucycles.c randombytes.c gaussian_ref.c
HEADERS256 = config.h params.h poly256.h ntt256.h reduce256.h symmetric.h randombytes.h fastrandombytes.h test/cpucycles.h gaussian_ref.h
KECCAK_SOURCES256 = $(SOURCES256) fips202.c symmetric-shake.c
KECCAK_HEADERS256 = $(HEADERS256) fips202.h
BENCH_HEADERS = test/bench.h test/cpu_freq_runtime.h
BENCH_SOURCES = test/bench.c test/cpu_freq_runtime.c

BFLAGS = -O3 -fomit-frame-pointer -march=native -flto -ftree-vectorize $(DEFINES)

ifdef DKG
    BFLAGS += -DBENCH_DKG
	BENCH_TARGET = test/dkg_$(USERS)_$(THRESHOLD)
endif

ifdef PREENC
    BFLAGS += -DBENCH_AS_PREENC
	BENCH_TARGET = test/preenc_$(USERS)_$(THRESHOLD)
endif

ifdef ENCRYPT
    BFLAGS += -DBENCH_AS_ENCRYPT
	BENCH_TARGET = test/encrypt_$(USERS)_$(THRESHOLD)
endif

ifdef AS2
    BFLAGS += -DBENCH_AS_2
	BENCH_TARGET = test/as2_$(USERS)_$(THRESHOLD)
endif


UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	GMP_DIR := /opt/homebrew/opt/gmp
	CFLAGS += -I$(GMP_DIR)/include
	LDFLAGS += -L$(GMP_DIR)/lib
endif

.PHONY: all speed shared clean

ghkss: \
  test/ghkss256_$(USERS)_$(THRESHOLD)

bench:\
  $(BENCH_TARGET)

test/ghkss256_$(USERS)_$(THRESHOLD): test/test_GHKSS256.c GHKSS256.c GHKSS256.h $(KECCAK_SOURCES256) $(KECCAK_HEADERS256)
	$(CC) $(CFLAGS) \
	  -o $@ test/test_GHKSS256.c GHKSS256.c $(KECCAK_SOURCES256) -lgmp -msse4.1 -maes -lm

$(BENCH_TARGET): test/benchmark_ghkss.c GHKSS256.c GHKSS256.h $(KECCAK_SOURCES256) $(KECCAK_HEADERS256) $(BENCH_HEADERS) $(BENCH_SOURCES)
	$(CC) $(BFLAGS) $(LDFLAGS) \
	  -o $@ test/benchmark_ghkss.c GHKSS256.c $(KECCAK_SOURCES256) $(BENCH_SOURCES) -lgmp -lm




clean:
	rm -f *~ test/*~ *.gcno *.gcda *.lcov
	rm -f libpqcrystals_dilithium2_ref.so
	rm -f libpqcrystals_dilithium3_ref.so
	rm -f libpqcrystals_dilithium5_ref.so
	rm -f libpqcrystals_fips202_ref.so
	rm -f test/test_dilithium2
	rm -f test/test_dilithium3
	rm -f test/test_dilithium5
	rm -f test/test_vectors2
	rm -f test/test_vectors3
	rm -f test/test_vectors5
	rm -f test/test_speed2
	rm -f test/test_speed3
	rm -f test/test_speed5
	rm -f test/test_mul
	rm -f nistkat/PQCgenKAT_sign2
	rm -f nistkat/PQCgenKAT_sign3
	rm -f nistkat/PQCgenKAT_sign5
	rm -f test/ghkss
	rm -f test/ghkss256*
	rm -f test/gaussian
	rm -f test/speed
	rm -f test/dkg*
	rm -f test/preenc*
	rm -f test/encrypt*
	rm -f test/as2*

