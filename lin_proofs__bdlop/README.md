# Olingo lin proofs

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