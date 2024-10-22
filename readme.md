## MPC-in-the-Head Without Repetition

This repository is the implementation of the MPC-in-the-Head-based NIZKAoK for the (modular)LWE instance. Original Paper: [IEEE](https://www.computer.org/csdl/proceedings-article/sp/2024/313000a157/1Ub24hrQXle) or [eprint](https://eprint.iacr.org/2024/1591)

- the folder `/kyber` contains the Kyber reference implementation from [Kyber Code](https://github.com/pq-crystals/kyber).
- the folder `/ntl` contains the NTL library from [NTL](https://www.shoup.net/ntl/) where you need to complie the static or dynamic linked library of NTL and set the path in `CMakeList.txt`before running the code.
- the `kosk.hpp` and `kosk.cpp` are the core APIs for the knowledge of secret key (KOSK) protocol for kyber KEM.
  - `kyber_verifiable_keygen` function returns a pair of public key and secret key (can be directly used for Kyber encapsulation and decapsulation api, see the `crypto_kem_enc` and `crypto_kem_dec` in `/kyber/kem` ) and proof (encoded into an array in uint_8) of its underlying LWE relation.
  - `kyber_kosk_verify` function verifies the proof of the underlying LWE relation.
- This implementation is executed in single-threaded mode. For multi-threaded mode, the experiment in paper using the multi-thread computation from NTL library.

> Note:
> The original version of code is stored in the `old` branch which has been used to run benchmarks in our papers.
> The code in the `main` branch is the updated version for better performance.
