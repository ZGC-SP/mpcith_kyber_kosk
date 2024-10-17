
#include "kosk.hpp"

void kyber_keygen(kyber_keypair *keypair,
                  mlwe_inst *raw_key) {
    uint8_t buf[2 * KYBER_SYMBYTES];
    const uint8_t * public_seed = buf;
    const uint8_t * noise_seed = buf + KYBER_SYMBYTES;
    uint8_t nonce = 0;
    polyvec a[KYBER_K], e, pk_raw, sk_raw;

    randombytes(buf, 2 * KYBER_SYMBYTES);
    buf[KYBER_SYMBYTES] = KYBER_K;
    sha3_512(buf, buf, KYBER_SYMBYTES + 1); // hash_g

    gen_matrix(a, public_seed, 0);
    for(size_t i = 0;i < KYBER_K; i++)
        poly_getnoise_eta1(&sk_raw.vec[i], noise_seed, nonce++);
    for(size_t i = 0;i < KYBER_K; i++)
        poly_getnoise_eta1(&e.vec[i], noise_seed, nonce++);

    // set raw key material
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_K; ++j) {
            for (size_t k = 0; k < KYBER_N; ++k) {
                raw_key->A[i].vec[j].coeffs[k] = a[i].vec[j].coeffs[k];
            }
        }
    }

    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_N; ++j) {
            raw_key->s.vec[i].coeffs[j] = sk_raw.vec[i].coeffs[j];
            raw_key->e.vec[i].coeffs[j] = e.vec[i].coeffs[j];
        }
    }

    // continue key generation
    polyvec_ntt(&sk_raw);
    polyvec_ntt(&e);

    for (size_t i = 0; i < KYBER_K; ++i) {
        polyvec_basemul_acc_montgomery(&pk_raw.vec[i], &a[i], &sk_raw);
        poly_tomont(&pk_raw.vec[i]);
    }

    polyvec_add(&pk_raw, &pk_raw, &e);
    polyvec_reduce(&pk_raw);

    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_N; ++j) {
            raw_key->t.vec[i].coeffs[j] = pk_raw.vec[i].coeffs[j];
        }
    }

    // pack_pk(keypair->pk, &pk_raw);
    polyvec_tobytes(keypair->pk, &pk_raw);
    memcpy(keypair->pk + KYBER_POLYVECBYTES, public_seed, KYBER_SYMBYTES);

    // sk_kem = (sk_pke || pk || H(pk) || z)
    // pack_sk(keypair->sk, &sk_raw);
    polyvec_tobytes(keypair->sk, &sk_raw);
    memcpy(keypair->sk + KYBER_INDCPA_SECRETKEYBYTES, keypair->pk, KYBER_PUBLICKEYBYTES);
    sha3_256(keypair->sk + KYBER_SECRETKEYBYTES - (2 * KYBER_SYMBYTES),
             keypair->pk,
             KYBER_PUBLICKEYBYTES);
    memcpy(keypair->sk + KYBER_SECRETKEYBYTES - KYBER_SYMBYTES,
           buf + KYBER_SYMBYTES,
           KYBER_SYMBYTES);
}

void kyber_verifiable_keygen(kyber_keypair *keypair,
                             uint8_t *pi) {
    mlwe_inst raw_mlwe;
    kyber_keygen(keypair, &raw_mlwe);

    mpcith_randomness rand0;
    mpcith_range_proof rand1;

    prepare_randomness(&rand0);
    prepare_range_proof(&rand1);

    mpcith_proof proof;
    prove(&proof, &raw_mlwe, &rand0, &rand1);
    encode_mpcith_proof(pi, &proof);
}

bool kyber_kosk_verify(const uint8_t *pi,
                       const uint8_t *pk) {
    // decode pk
    polyvec A[KYBER_K], t;
    uint8_t seed[KYBER_SYMBYTES];

    polyvec_frombytes(&t, pk);
    memcpy(seed, pk + KYBER_POLYVECBYTES, KYBER_SYMBYTES);

    gen_matrix(A, seed, 0);

    mlwe_inst raw_mlwe;
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_K; ++j) {
            for (size_t k = 0; k < KYBER_N; ++k) {
                raw_mlwe.A[i].vec[j].coeffs[k] = A[i].vec[j].coeffs[k];
            }
        }
    }
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_N; ++j) {
            raw_mlwe.t.vec[i].coeffs[j] = t.vec[i].coeffs[j];
        }
    }

    mpcith_proof proof;
    decode_mpcith_proof(&proof, pi);

    return verify(&proof, &raw_mlwe);
}