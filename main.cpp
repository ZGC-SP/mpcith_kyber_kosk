#include <iostream>

extern "C" {
#include "kyber/kem.h"
#include "kyber/indcpa.h"
}
#include "mlwe_prover.hpp"
#include "kosk.hpp"
#include "params.hpp"

using namespace NTL;
int main() {
    /**
     * MLWE prover.
     */
    printf("=== mlwe prover test === \n");
    // test share ranges
    clock_t start_pre, end_pre;
    start_pre = clock();

    mpcith_randomness rand;
    mpcith_range_proof eta_shares;
    prepare_randomness(&rand);
    prepare_range_proof(&eta_shares);

    end_pre = clock();
    double cpu_time_used_pre = ((double) (end_pre - start_pre)) / CLOCKS_PER_SEC;
    printf(">>> preprocessing time used: %f s\n", cpu_time_used_pre);

    mlwe_inst raw_sec; kyber_keypair kp;
    kyber_keygen(&kp, &raw_sec);

    clock_t start_proof, end_proof;
    start_proof = clock();

    mpcith_proof pi;
    prove(&pi, &raw_sec, &rand, &eta_shares);

    end_proof = clock();
    double cpu_time_used_proof = ((double) (end_proof - start_proof)) / CLOCKS_PER_SEC;
    printf(">>> prove time used: %f s\n", cpu_time_used_proof);
    printf(">>> total prove time used: %f s\n", cpu_time_used_pre + cpu_time_used_proof);

    // verify phase:
    clock_t start_verify = clock();

    bool res = verify(&pi, &raw_sec);

    clock_t end_verify = clock();
    double cpu_time_used_verify = ((double) (end_verify - start_verify)) / CLOCKS_PER_SEC;
    printf(">>> verify time used: %f s \n", cpu_time_used_verify);

    if (res) {
        printf("[result] mlwe verify success\n");
    } else {
        printf("[result] mlwe verify failed\n");
    }

    printf("[proof size] %lu kilobytes\n", MPCITH_PROOF_SIZE / 1024);

    printf("\n");

    /**
     * Kyber verifiable key generation and Knowldege of Secret Key (KOSK) verifier.
     */

    printf("=== kyber verifiable keygen & KOSK === \n");
    clock_t start_kosk_keygen, end_kosk_keygen;

    kyber_keypair keypair;
    uint8_t kosk_pi[MPCITH_PROOF_SIZE] = {0};

    start_kosk_keygen = clock();
    kyber_verifiable_keygen(&keypair, kosk_pi);
    end_kosk_keygen = clock();

    double cpu_time_used_kosk_keygen = ((double) (end_kosk_keygen - start_kosk_keygen)) / CLOCKS_PER_SEC;
    printf(">>> kyber verifiable keygen (keygen + preprocess + prove) time used: %f s\n", cpu_time_used_kosk_keygen);

    // verify
    clock_t start_kosk_verify, end_kosk_verify;

    start_kosk_verify = clock();
    bool res2 = kyber_kosk_verify(kosk_pi, keypair.pk);
    end_kosk_verify = clock();

    if (res2) {
        printf("[result] kyber kosk verify success\n");
    } else {
        printf("[result] kyber kosk verify failed\n");
    }

    double cpu_time_used_kosk_verify = ((double) (end_kosk_verify - start_kosk_verify)) / CLOCKS_PER_SEC;
    printf(">>> kyber kosk verify time used: %f s\n", cpu_time_used_kosk_verify);
    printf("\n");

    // test encap and decap
    printf("=== kyber kem: encaps & decaps === \n");
    uint8_t ss2[KYBER_SSBYTES], ss_dec2[KYBER_SSBYTES], ct2[KYBER_CIPHERTEXTBYTES];

    crypto_kem_enc(ct2, ss2, keypair.pk);
    crypto_kem_dec(ss_dec2, ct2, keypair.sk);
    // check ss2 and ss_dec2
    for (int i = 0; i < KYBER_SSBYTES; ++i) {
        if (ss2[i] != ss_dec2[i]) {
            printf("\t decapsulated ss is not the same as the encapsulated\n");
            break;
        }
        if (i == KYBER_SSBYTES - 1) {
            printf("[result] decapsulated ss is the same as the encapsulated\n");
        }
    }


    return 0;
}
