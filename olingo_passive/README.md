# Passive Olingo
The code for the passively secure Olingo scheme is based on the Dilithium code base.
For the implementation we have used Dilithium's polynomials and ntt functions with GMP integers instead of the int64_t type in our implementation.

Passive signing benchmarks can be generated with `make ghkss`, from `test_GHKSS256.c`.
Passive DKG and keygen benchmarks can be generated with `make bench`, from `benchmark_ghkss.c`.

```
make bench USERS=<2-1024> THRESHOLD=<1-1023> DKG=1        # Builds DKG benchmarking
make bench USERS=<2-1024> THRESHOLD=<1-1023> PREENC=1     # Builds first part of KGenS
make bench USERS=<2-1024> THRESHOLD=<1-1023> ENCRYPT=1    # Builds second part of KGenS
make bench USERS=<2-1024> THRESHOLD=<1-1023> AS2=1        # Builds third part of KGenS
make ghkss USERS=<2-1024> THRESHOLD=<1-1023>              # Builds signing benchmarks
```

To build with specific threshold, build using f.ex.:
`make bench DKG USERS=4 THRESHOLD=2`

Or, you can build all benchmarking files for one set of users and threshold using the bash script `build.sh`:
```
bash build.sh --users=<int> --threshold=<int>
```


## Build instructions

The implementations contain several test and benchmarking programs and a Makefile to facilitate compilation.

### Prerequisites

Some of the test programs require [OpenSSL](https://openssl.org). If the OpenSSL header files and/or shared libraries do not lie in one of the standard locations on your system, it is necessary to specify their location via compiler and linker flags in the environment variables `CFLAGS`, `NISTFLAGS`, and `LDFLAGS`.

For example, on macOS you can install OpenSSL via [Homebrew](https://brew.sh) by running
```sh
brew install openssl
```
Then, run
```sh
export CFLAGS="-I/opt/homebrew/opt/openssl@1.1/include"
export NISTFLAGS="-I/opt/homebrew/opt/openssl@1.1/include"
export LDFLAGS="-L/opt/homebrew/opt/openssl@1.1/lib"
```
before compilation to add the OpenSSL header and library locations to the respective search paths.

### This code is based on a skeleton of the Dilithium C implementation (https://github.com/pq-crystals/dilithium) using their reference implementation.