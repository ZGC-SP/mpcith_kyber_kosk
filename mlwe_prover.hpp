
#ifndef MPCITH_KYBER_KOSK_MLWE_PROVER_HPP
#define MPCITH_KYBER_KOSK_MLWE_PROVER_HPP

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
using namespace NTL;

#include "ss.hpp"
#include "params.hpp"

extern "C" {
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#include "kyber/symmetric.h"
#include "kyber/params.h"
#include "kyber/polyvec.h"
#include "kyber/indcpa.h"
#include "utils//gf3329.h"
#include "utils/precomputed_kyber.h"
}

// constant for mpc-in-the-head proofs
#define MPCITH_PROOF_SIZE sizeof(mpcith_proof)
#define MPCITH_PRE_RANDOMNESS_SIZE sizeof(mpcith_randomness) + sizeof(mpcith_range_proof)
//

typedef struct {
    polyvec A[KYBER_K], t;
    polyvec s, e;
} mlwe_inst;

typedef struct {
    uint16_t f[MPCITH_K + MPCITH_V + 1][KYBER_N];
    uint16_t NTT_f[MPCITH_K + MPCITH_V + 1][KYBER_N];
    share_vec f_shares[MPCITH_K + MPCITH_V + 1];
    share_vec NTT_f_shares[MPCITH_K + MPCITH_V + 1];
} mpcith_randomness;

typedef struct {
    share_vec s_eta_shares[KYBER_K][KYBER_ETA1 * 2 + 1];
    share_vec e_eta_shares[KYBER_K][KYBER_ETA1 * 2 + 1];
} mpcith_range_proof;
void encode_preprocessed_randomness(uint8_t *buf,
                                    const mpcith_randomness *rand,
                                    const mpcith_range_proof *eta_shares);
void decode_preprocessed_randomness(mpcith_randomness *rand,
                                    mpcith_range_proof *eta_shares,
                                    const uint8_t *buf);

typedef struct {
    // preprocessing:
    uint16_t f_shares[MPCITH_T][MPCITH_K + MPCITH_V + 1], NTT_f_shares[MPCITH_T][MPCITH_K + MPCITH_V + 1];
    uint16_t beta_shares[MPCITH_N - MPCITH_T][MPCITH_K], gamma_shares[MPCITH_N - MPCITH_T][MPCITH_K];
    uint8_t Tcomm[MPCITH_N - MPCITH_T][KYBER_SYMBYTES];
    // online:
    // (1) Relation
    uint16_t I[MPCITH_T];
    uint16_t s_shares[MPCITH_T][KYBER_K], e_shares[MPCITH_T][KYBER_K], t_shares[MPCITH_N - MPCITH_T][KYBER_K];
    uint16_t NTT_s_shares[MPCITH_T][KYBER_K], NTT_e_shares[MPCITH_T][KYBER_K];
    uint16_t NTT_Ar_shares[MPCITH_T][KYBER_K], NTT_As_shares[MPCITH_T][KYBER_K];
    uint16_t sr_shares[MPCITH_N - MPCITH_T][KYBER_K], er_shares[MPCITH_N - MPCITH_T][KYBER_K];
    // (2) short proofs
    uint16_t s_eta_shares[MPCITH_N - MPCITH_T][KYBER_K][KYBER_ETA1 * 2 + 1], e_eta_shares[MPCITH_N - MPCITH_T][KYBER_K][KYBER_ETA1 * 2 + 1];
    uint16_t s_sub_eta_shares[MPCITH_T][KYBER_K][KYBER_ETA1 * 2 + 1], e_sub_eta_shares[MPCITH_T][KYBER_K][KYBER_ETA1 * 2 + 1];
    uint16_t z_s_ddeg_shares[MPCITH_T][KYBER_K][KYBER_ETA1 * 2], z_e_ddeg_shares[MPCITH_T][KYBER_K][KYBER_ETA1 * 2];
    uint16_t u_s_2ddeg_shares[MPCITH_N - MPCITH_T][KYBER_K][KYBER_ETA1 * 2], u_e_2ddeg_shares[MPCITH_N - MPCITH_T][KYBER_K][KYBER_ETA1 * 2];
    uint8_t comm[MPCITH_N - MPCITH_T][KYBER_SYMBYTES];
} mpcith_proof;

void encode_mpcith_proof(uint8_t *buf,
                         const mpcith_proof *pi);
void decode_mpcith_proof(mpcith_proof *pi,
                         const uint8_t *buf);

void prepare_randomness(mpcith_randomness *rand);
void prepare_range_proof(mpcith_range_proof *eta_shares);

typedef struct {
    uint8_t comm[KYBER_SYMBYTES];
    uint16_t s_sh[KYBER_K], e_sh[KYBER_K];
    uint16_t f_sh[MPCITH_K + MPCITH_V + 1], Tf_sh[MPCITH_K + MPCITH_V + 1];
    uint16_t beta[MPCITH_K], gamma[MPCITH_K];
    // states in gates proof
    uint16_t sr_sh[KYBER_K], er_sh[KYBER_K];
    uint16_t s_ddeg_sh[KYBER_K][2 * KYBER_ETA1], e_ddeg_sh[KYBER_K][2 * KYBER_ETA1];
    uint16_t s_zero_sh[KYBER_K][2 * KYBER_ETA1], e_zero_sh[KYBER_K][2 * KYBER_ETA1];
} mpcith_vp_state;

void prove(mpcith_proof *pi,
           const mlwe_inst *mlwe,
           const mpcith_randomness *rand,
           const mpcith_range_proof *eta_share);




#endif //MPCITH_KYBER_KOSK_MLWE_PROVER_HPP
