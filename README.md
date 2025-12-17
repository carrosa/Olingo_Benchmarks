# Benchmarking Olingo

Passive signing benchmarks can be generated with `make ghkss`, from `test_GHKSS256.c`.
Passive DKG and keygen benchmarks can be generated with `make bench`, from `benchmark_ghkss.c`.

To benchmark with different thresholds and total parties, one needs to change the THRESHOLD and USERS variables in params.h. Keep in mind benchmarking DKG with large THRESHOLD and USERS is slow and requires quite a bit of memory.

We simulate 1 party only, meaning we use dummy data for the benchmarks where applicable to simulate having multiple parties.

This README file is an aggregation of the README's for each of the subfolders in this repository. 

## Passive Olingo
We build upon dilithiums codebase using GMP for larger integers.
The code for the passively secure Olingo scheme is based on the Dilithium code base.
For the implementation we have used Dilithium's polynomials and ntt functions with GMP integers instead of the int64_t type in our implementation.

Passive signing benchmarks can be generated with `make ghkss`, from `test_GHKSS256.c`.
Passive DKG and keygen benchmarks can be generated with `make bench`, from `benchmark_ghkss.c`.
See below or in the README.md in the folder `olingo_passive` for how to build.

**Notice:** The Makefile and build script build.sh is found inside the `ref` folder inside the `olingo_passive` folder.

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

### NTT
To compute zeta values for NTT, we used the script `compute_zetas.py` which has dependencies:
- sage
- argparse
- numpy

The information about requirements to run is found inside the folder `olingo_passive` or below.

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

## LIN proofs:

We adapt https://github.com/dfaranha/lattice-verifiable-mixnet for BDLOP commitments and linear proofs. See its readme for more info on how to build.
The original code for BDLOP and linearity proofs are based on the Pi_LIN BDLOP'18 paper, which we also use for our LIN proofs.

Information about requirements to run is found inside the folder `lin_proofs__bdlop` or below.

### Olingo lin proofs

Depedencies are the [NFLlib](https://github.com/quarkslab/NFLlib) and [FLINT](https://flintlib.org/doc/) 2.9 libraries.
NFLLib is already included in this repository, but instructions for installing its dependencies can be found in the link above.
FLINT is usually included in package managers and can be easily installed in most systems out there.

### Building dependencies

To build NFLLib, run the following inside a cloned version of this repository:

```
$ mkdir deps
$ cd deps
$ cmake ../NFLlib -DCMAKE_BUILD_TYPE=Release -DNFL_OPTIMIZED=ON
$ make
$ make test
```
### Building and running the code

For building the Olingo LIN proof benchmarks, run `make olingo` inside the source directory. This will build the binary `olingo` to benchmark the LIN proofs in our paper.

__WARNING__: This is an academic proof of concept, and in particular has not received code review. This implementation is NOT ready for any type of production use.


### This code is originally forked from the code accompannying the paper "Verifiable Mix-Nets and Distributed Decryption for Voting from Lattice-Based Assumptions".
[mixnets](https://github.com/dfaranha/lattice-verifiable-mixnet)

## Olingo lazer proof benchmarks:

We use the python interface for Lazer for LNP style proofs.
The folders for our benchmarks:

* python/bnd_E_prf - norm bound proof of E with sigma_tdec (norm bound proof of pi_dsi)
* python/dkg_prf1 - prove norm bound of S for dkg (norm bound proof for pi_{KGen_E}
* python/pi_r - proving pi_r (sigma=2^38, otherwise same as pi_si)
* python/pi_si - proving pi_si (sigma=2^20, otherwise same as pi_r)

To build each of the benchmarks, run `make` in their respective folders.

More information about requirements to run is found inside the folder `olingo__lazer`, in the file `lazer_getting_started.html` and in the github repository for lazer: https://github.com/lazer-crypto/lazer.


# Estimating Olingo Security
The scripts used to estimate hardness with our selected parameters can be found in the folder `lattice-estimator`.
They depend on the lattice-estimator [github.com/malb/lattice-estimator](https://github.com/malb/lattice-estimator).

# A warning to developers:
__WARNING__: This is an academic proof of concept/benchmarking, and in particular has not received code review. This implementation is NOT ready for any type of production use.
