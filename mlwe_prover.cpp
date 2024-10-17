
#include "mlwe_prover.hpp"

void prepare_randomness(mpcith_randomness *rand) {
    // generate k + v + 1 random vectors
    uint8_t seeds[MPCITH_K + MPCITH_V + 1][KYBER_SYMBYTES];
    uint8_t prand_bytes[MPCITH_K + MPCITH_V + 1][2 * KYBER_N]; // 2 bytes for elements over GF(3329)
    for (size_t i = 0; i < MPCITH_K + MPCITH_V + 1; ++i) {
        randombytes(seeds[i], KYBER_SYMBYTES);
        kyber_shake256_prf(prand_bytes[i], 2 * KYBER_N, seeds[i], i);
        for (size_t j = 0; j < KYBER_N; ++j) {
            rand->f[i][j] = ((prand_bytes[i][2 * j] << 8) | prand_bytes[i][2 * j + 1]) % KYBER_Q;
        }
    }

    // evaluate NTT(f)
    poly NTT_f[MPCITH_K + MPCITH_V + 1];
    for (size_t i = 0; i < MPCITH_K + MPCITH_V + 1; ++i) {
        for (size_t j = 0; j < KYBER_N; ++j) {
            NTT_f[i].coeffs[j] = decode_from_gf3329(rand->f[i][j]);
        }
        poly_ntt(&NTT_f[i]);
        for (size_t j = 0; j < KYBER_N; ++j) {
            rand->NTT_f[i][j] = encode_to_gf3329(NTT_f[i].coeffs[j]);
        }
    }

    // share the randomness
    for (size_t i = 0; i < MPCITH_K + MPCITH_V + 1; ++i) {
        // share f
        secret_vec f_sec;
        memcpy(f_sec.secret, rand->f[i], KYBER_N * sizeof(uint16_t));
        share_secrets_ddeg(&rand->f_shares[i], &f_sec);
        // share NTT(f)
        secret_vec NTT_f_sec;
        memcpy(NTT_f_sec.secret, rand->NTT_f[i], KYBER_N * sizeof(uint16_t));
        share_secrets_ddeg(&rand->NTT_f_shares[i], &NTT_f_sec);
    }
}

void prepare_range_proof(mpcith_range_proof *eta_share) {
    uint16_t eta_sec[2 * KYBER_ETA1 + 1][KYBER_N];
    for (int16_t i = (-1) * KYBER_ETA1; i <= KYBER_ETA1; ++i) {
        uint16_t e = encode_to_gf3329(i);
        for (size_t j = 0; j < KYBER_N; ++j) {
            eta_sec[i + KYBER_ETA1][j] = e;
        }
    }

    // share eta range
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < 2 * KYBER_ETA1 + 1; ++j) {
            secret_vec eta_sec_vec;
            memcpy(eta_sec_vec.secret, eta_sec[j], KYBER_N * sizeof(uint16_t));
            share_secrets_ddeg(&eta_share->s_eta_shares[i][j], &eta_sec_vec);
            share_secrets_ddeg(&eta_share->e_eta_shares[i][j], &eta_sec_vec);
        }
    }
}

void encode_preprocessed_randomness(uint8_t *buf,
                                    const mpcith_randomness *rand,
                                    const mpcith_range_proof *eta_shares) {
    memcpy(buf, rand, sizeof(mpcith_randomness));
    memcpy(buf + sizeof(mpcith_randomness), eta_shares, sizeof(mpcith_range_proof));

}
void decode_preprocessed_randomness(mpcith_randomness *rand,
                                    mpcith_range_proof *eta_shares,
                                    const uint8_t *buf) {
    for (size_t i = 0; i < MPCITH_K + MPCITH_V + 1; ++i) {
        memcpy(rand->f[i] , ((mpcith_randomness *) buf)->f[i], KYBER_N * sizeof(uint16_t));
        memcpy(rand->NTT_f[i], ((mpcith_randomness *) buf)->NTT_f[i], KYBER_N * sizeof(uint16_t));
        memcpy(rand->f_shares[i].share_x, ((mpcith_randomness *) buf)->f_shares[i].share_x, MPCITH_N * sizeof(uint16_t));
        memcpy(rand->f_shares[i].share_y, ((mpcith_randomness *) buf)->f_shares[i].share_y, MPCITH_N * sizeof(uint16_t));
        memcpy(rand->NTT_f_shares[i].share_x, ((mpcith_randomness *) buf)->NTT_f_shares[i].share_x, MPCITH_N * sizeof(uint16_t));
        memcpy(rand->NTT_f_shares[i].share_y, ((mpcith_randomness *) buf)->NTT_f_shares[i].share_y, MPCITH_N * sizeof(uint16_t));
    }
}

void prove(mpcith_proof *pi,
           const mlwe_inst *mlwe,
           const mpcith_randomness *rand,
           const mpcith_range_proof *eta_share) {
    /**
     * preprocessing: prepare the r and NTT(r) shared pairs
     */
    // encode s, e in mlwe over GF(3329)
    secret_vec s_sec[KYBER_K], e_sec[KYBER_K];
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_N; ++j) {
            s_sec[i].secret[j] = encode_to_gf3329(mlwe->s.vec[i].coeffs[j]);
            e_sec[i].secret[j] = encode_to_gf3329(mlwe->e.vec[i].coeffs[j]);
        }
    }
    // share s, e
    share_vec s_shares[KYBER_K], e_shares[KYBER_K];
    for (size_t i = 0; i < KYBER_K; ++i) {
        share_secrets_ddeg(&s_shares[i], &s_sec[i]);
        share_secrets_ddeg(&e_shares[i], &e_sec[i]);
    }

    mpcith_vp_state V[MPCITH_N];
    for (size_t i = 0; i < MPCITH_N; ++i) {
        for (size_t j = 0; j < KYBER_K; ++j) {
            V[i].s_sh[j] = s_shares[j].share_y[i];
            V[i].e_sh[j] = e_shares[j].share_y[i];
        }
        for (size_t j = 0; j < MPCITH_K + MPCITH_V + 1; ++j) {
            V[i].f_sh[j] = rand->f_shares[j].share_y[i];
            V[i].Tf_sh[j] = rand->NTT_f_shares[j].share_y[i];
        }
    }

    // commit: T_com_i = Hash(s_sh_i, e_sh_i, f_sh_i, T_sh_i)
    for (size_t i = 0; i < MPCITH_N; ++i) {
        uint16_t ctxt[(KYBER_K + MPCITH_K + MPCITH_V + 1) * 2];
        for (size_t j = 0; j < KYBER_K; ++j) {
            ctxt[j] = V[i].s_sh[j];
            ctxt[j + KYBER_K] = V[i].e_sh[j];
        }
        for (size_t j = 0; j < MPCITH_K + MPCITH_V + 1; ++j) {
            ctxt[j + (KYBER_K * 2)] = V[i].f_sh[j];
            ctxt[j + (KYBER_K * 2) + (MPCITH_K + MPCITH_V + 1)] = V[i].Tf_sh[j];
        }
        sha3_256(V[i].comm, (const uint8_t *) ctxt, sizeof(ctxt));
    }

    // prepare alpha list for r
    uint8_t alpha_[sizeof(uint16_t) * (MPCITH_K + MPCITH_V)];
    uint8_t Tcomm[MPCITH_N * KYBER_SYMBYTES], digest[KYBER_SYMBYTES];
    for (size_t i = 0; i < MPCITH_N; ++i) {
        memcpy(Tcomm + i * KYBER_SYMBYTES, V[i].comm, KYBER_SYMBYTES);
    }
    sha3_256(digest, Tcomm, MPCITH_N * KYBER_SYMBYTES);
    kyber_shake256_prf(alpha_, sizeof(uint16_t) * (MPCITH_K + MPCITH_V), digest, 1);

    uint16_t alpha[MPCITH_K + MPCITH_V];
    // encode uint8 to GF(3329)
    for (size_t i = 0; i < MPCITH_K + MPCITH_V; ++i) {
        alpha[i] = ((alpha_[2 * i] << 8) | alpha_[2 * i + 1]) % KYBER_Q;
    }

    uint16_t pow_table[MPCITH_K + MPCITH_V][MPCITH_K + MPCITH_V + 1];
    for (size_t i = 0; i < MPCITH_K + MPCITH_V; ++i) {
        pow_table[i][0] = 1;
        pow_table[i][1] = alpha[i];
    }
    for (size_t i = 0; i < MPCITH_K + MPCITH_V; ++i) {
        for (size_t j = 2; j < MPCITH_K + MPCITH_V + 1; ++j) {
            pow_table[i][j] = gf3329_mul(pow_table[i][j - 1], alpha[i]);
        }
    }

    // evaluate:
    // beta(i) = \sum_{j=0}^{k+v} alpha^j * f_sh_j(i)
    // gamma(i) = \sum_{j=0}^{k+v} alpha^j * Tf_sh_j(i)

    for (size_t i = 0; i < MPCITH_N; ++i) {
        for (size_t j = 0; j < MPCITH_K; ++j) {
            V[i].beta[j] = 0;
            for (size_t k = 0; k < MPCITH_K + MPCITH_V + 1; ++k) {
                if (k == 0) {
                    V[i].beta[j] = V[i].f_sh[k];
                } else {
                    V[i].beta[j] = gf3329_add(V[i].beta[j], gf3329_mul(pow_table[j][k], V[i].f_sh[k]));
                }
            }
            V[i].gamma[j] = 0;
            for (size_t k = 0; k < MPCITH_K + MPCITH_V + 1; ++k) {
                if (k == 0) {
                    V[i].gamma[j] = V[i].Tf_sh[k];
                } else {
                    V[i].gamma[j] = gf3329_add(V[i].gamma[j], gf3329_mul(pow_table[j][k], V[i].Tf_sh[k]));
                }
            }
        }
    }

    // evaluate r, NTT_r
    uint16_t r[MPCITH_N][MPCITH_V], NTT_r[MPCITH_N][MPCITH_V];
    for (size_t i = 0; i < MPCITH_N; ++i) {
        for (size_t j = MPCITH_K; j < MPCITH_V + MPCITH_K; ++j) {
            r[i][j - MPCITH_K] = 0;
            for (size_t k = 0; k < MPCITH_K + MPCITH_V + 1; ++k) {
                if (k == 0) {
                    r[i][j - MPCITH_K] = V[i].f_sh[k + MPCITH_K + 1];
                } else {
                    r[i][j - MPCITH_K] = gf3329_add(r[i][j - MPCITH_K],
                                                    gf3329_mul(pow_table[j][k], V[i].f_sh[k]));
                }
            }
            NTT_r[i][j - MPCITH_K] = 0;
            for (size_t k = 0; k < MPCITH_K + MPCITH_V + 1; ++k) {
                if (k == 0) {
                    NTT_r[i][j - MPCITH_K] = V[i].Tf_sh[k + MPCITH_K + 1];
                } else {
                    NTT_r[i][j - MPCITH_K] = gf3329_add(NTT_r[i][j - MPCITH_K],
                                                        gf3329_mul(pow_table[j][k], V[i].Tf_sh[k]));
                }
            }
        }
    }

    // encode to shared type
    share_vec r_share[MPCITH_V], NTT_r_share[MPCITH_V];
    for (size_t i = 0; i < MPCITH_V; ++i) {
        for (size_t j = 0; j < MPCITH_N; ++j) {
            r_share[i].share_x[j] = j + KYBER_N;
            r_share[i].share_y[j] = r[j][i];
            NTT_r_share[i].share_x[j] = j + KYBER_N;
            NTT_r_share[i].share_y[j] = NTT_r[j][i];
        }
    }

    /**
     * A * s + e: relation proof: A * NTT(s) + NTT(e) = NTT(t)
     */

    // add: [s] + [r1], [e] + [r2]
    // recons: s' = s + r1, e' = e + r2
    share_vec sr_shares[KYBER_K], er_shares[KYBER_K];
    secret_vec sr_rec[KYBER_K], er_rec[KYBER_K];
    uint16_t sr_rnd[KYBER_K][DEG_D + 1], er_rnd[KYBER_K][DEG_D + 1];
    for (size_t i = 0; i < KYBER_K; ++i) {
        // s:
        shares_add(&sr_shares[i], &s_shares[i], &r_share[i]);
        recon_secrets_ddeg(&sr_rec[i], &sr_shares[i]);

        // e:
        shares_add(&er_shares[i], &e_shares[i], &r_share[i + KYBER_K]);
        recon_secrets_ddeg(&er_rec[i], &er_shares[i]);

        for (size_t j = KYBER_N; j < DEG_D + 1; ++j) {
            sr_rnd[i][j] = sr_shares[i].share_y[j - KYBER_N];
            er_rnd[i][j] = er_shares[i].share_y[j - KYBER_N];
        }
    }
    // set shares in view
    for (size_t i = 0; i < MPCITH_N; ++i) {
        for (size_t j = 0; j < KYBER_K; ++j) {
            V[i].sr_sh[j] = sr_shares[j].share_y[i];
            V[i].er_sh[j] = er_shares[j].share_y[i];
        }
    }
    // recover the NTT(r1) for s
    secret_vec ntt_r1_sec[KYBER_K]; polyvec ntt_r1, s_cpy, e_cpy;
    for (size_t i = 0; i < KYBER_K; ++i) {
        recon_secrets_ddeg(&ntt_r1_sec[i], &r_share[i]);
        for (size_t j = 0; j < KYBER_N; ++j) {
            ntt_r1.vec[i].coeffs[j] = decode_from_gf3329(ntt_r1_sec[i].secret[j]);
            s_cpy.vec[i].coeffs[j] = mlwe->s.vec[i].coeffs[j];
            e_cpy.vec[i].coeffs[j] = mlwe->e.vec[i].coeffs[j];
        }
    }
    polyvec_ntt(&s_cpy);
    polyvec_ntt(&e_cpy);

    // do ntt over s' and e' compute: t' = A * s' + e'
    polyvec sr, er, ntt_Asr, ntt_Ar, ntt_As, ntt_t;
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_N; ++j) {
            sr.vec[i].coeffs[j] = decode_from_gf3329(sr_rec[i].secret[j]);
            er.vec[i].coeffs[j] = decode_from_gf3329(er_rec[i].secret[j]);
            // t (after ntt)
            ntt_t.vec[i].coeffs[j] = mlwe->t.vec[i].coeffs[j];
        }
    }
    polyvec_ntt(&sr);
    polyvec_ntt(&er);
    // set rnd
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_N; ++j) {
            sr_rnd[i][j] = encode_to_gf3329(sr.vec[i].coeffs[j]);
            er_rnd[i][j] = encode_to_gf3329(er.vec[i].coeffs[j]);
        }
    }

    for (size_t i = 0; i < KYBER_K; ++i) {
        // A * r
        polyvec_basemul_acc_montgomery(&ntt_Ar.vec[i], &mlwe->A[i], &ntt_r1);
        poly_tomont(&ntt_Ar.vec[i]);
        // A * s
        polyvec_basemul_acc_montgomery(&ntt_As.vec[i], &mlwe->A[i], &s_cpy);
        poly_tomont(&ntt_As.vec[i]);
        // A * (s + r)
        polyvec_basemul_acc_montgomery(&ntt_Asr.vec[i], &mlwe->A[i], &sr);
        poly_tomont(&ntt_Asr.vec[i]);
    }

    // re-share NTT(s + r) and NTT(e + r) with rnd.
    share_vec ntt_sr_shares[KYBER_K], ntt_er_shares[KYBER_K];
    share_vec ntt_s_shares[KYBER_K], ntt_e_shares[KYBER_K];
    share_vec ntt_Ar_shares[KYBER_K], ntt_As_shares[KYBER_K], ntt_Asr_shares[KYBER_K];
    share_vec ntt_t_shares[KYBER_K];
    for (size_t i = 0; i < KYBER_K; ++i) {
        // share: [NTT(s + r1)], [NTT(e + r2)] with original randomness
        recompute_share_secrets_ddeg(&ntt_sr_shares[i], sr_rnd[i]);
        recompute_share_secrets_ddeg(&ntt_er_shares[i], er_rnd[i]);

        // compute [NTT(s)] = [NTT(s + r1)] - [NTT(r1)]; [NTT(e)] = [NTT(e + r2)] - [NTT(r2)]
        shares_sub(&ntt_s_shares[i],  &ntt_sr_shares[i], &NTT_r_share[i]);
        shares_sub(&ntt_e_shares[i], &ntt_er_shares[i], &NTT_r_share[i + KYBER_K]);
    }
    // reshared As, Ar, A(s+r)
    uint16_t ntt_Asr_rnd[KYBER_K][DEG_D + 1]; secret_vec ntt_As_sec[KYBER_K];
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_N; ++j) {
            ntt_As_sec[i].secret[j] = encode_to_gf3329(ntt_As.vec[i].coeffs[j]);
            ntt_Asr_rnd[i][j] = encode_to_gf3329(ntt_Asr.vec[i].coeffs[j]);
        }
        for (size_t j = KYBER_N; j < DEG_D + 1; ++j) {
            ntt_Asr_rnd[i][j] = sr_rnd[i][j];
        }
        recompute_share_secrets_ddeg(&ntt_Asr_shares[i], ntt_Asr_rnd[i]);
        share_secrets_ddeg(&ntt_As_shares[i], &ntt_As_sec[i]);
        shares_sub(&ntt_Ar_shares[i], &ntt_Asr_shares[i], &ntt_As_shares[i]);
    }

    // add gate: A * s + e = t
    for (size_t i = 0; i < KYBER_K; ++i) {
        shares_add(&ntt_t_shares[i], &ntt_As_shares[i], &ntt_e_shares[i]);
    }

    // test recon NTT(t)
    secret_vec recon_t_sec;
    recon_secrets_ddeg(&recon_t_sec, &ntt_t_shares[0]);
    // encode to poly
    poly recon_t;
    for (size_t i = 0; i < KYBER_N; ++i) {
        recon_t.coeffs[i] = decode_from_gf3329(recon_t_sec.secret[i]);
    }

    /**
     * range proof: prove that s and e are in the range of [-eta, eta]
     */
    // s - eta, e - eta
    share_vec s_sub_eta_shares[KYBER_K][2 * KYBER_ETA1 + 1], e_sub_eta_shares[KYBER_K][2 * KYBER_ETA1 + 1];
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < 2 * KYBER_ETA1 + 1; ++j) {
            shares_sub(&s_sub_eta_shares[i][j], &s_shares[i], &eta_share->s_eta_shares[i][j]);
            shares_sub(&e_sub_eta_shares[i][j], &e_shares[i], &eta_share->e_eta_shares[i][j]);
        }
    }

    // evaluate mul gates for s:
    share_vec s_reduced_ddeg_shares[KYBER_K][2 * KYBER_ETA1], s_reduced_2ddeg_shares[KYBER_K][2 * KYBER_ETA1];
    share_vec e_reduced_ddeg_shares[KYBER_K][2 * KYBER_ETA1], e_reduced_2ddeg_shares[KYBER_K][2 * KYBER_ETA1];
    share_vec s_zero_2ddeg_shares[KYBER_K][2 * KYBER_ETA1], e_zero_2ddeg_shares[KYBER_K][2 * KYBER_ETA1];

    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < 2 * KYBER_ETA1; ++j) {
            secret_vec s_reduced_ddeg, e_reduced_ddeg;
            // mul:
            if (j == 0) {
                shares_mul(&s_reduced_2ddeg_shares[i][j],
                           &s_sub_eta_shares[i][j], &s_sub_eta_shares[i][j + 1]);
                shares_mul(&e_reduced_2ddeg_shares[i][j],
                           &e_sub_eta_shares[i][j], &e_sub_eta_shares[i][j + 1]);
            } else {
                shares_mul(&s_reduced_2ddeg_shares[i][j],
                           &s_reduced_ddeg_shares[i][j - 1], &s_sub_eta_shares[i][j + 1]);
                shares_mul(&e_reduced_2ddeg_shares[i][j],
                           &e_reduced_ddeg_shares[i][j - 1], &e_sub_eta_shares[i][j + 1]);
            }

            // reduce degree
            recon_secrets_2ddeg(&s_reduced_ddeg, &s_reduced_2ddeg_shares[i][j]);
            share_secrets_ddeg(&s_reduced_ddeg_shares[i][j], &s_reduced_ddeg);
            recon_secrets_2ddeg(&e_reduced_ddeg, &e_reduced_2ddeg_shares[i][j]);
            share_secrets_ddeg(&e_reduced_ddeg_shares[i][j], &e_reduced_ddeg);
        }
    }
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < 2 * KYBER_ETA1; ++j) {
            shares_sub(&s_zero_2ddeg_shares[i][j],
                       &s_reduced_2ddeg_shares[i][j], &s_reduced_ddeg_shares[i][j]);
            shares_sub(&e_zero_2ddeg_shares[i][j],
                       &e_reduced_2ddeg_shares[i][j], &e_reduced_ddeg_shares[i][j]);
        }
    }
    // set mul-gate shares in view
    for (size_t i = 0; i < MPCITH_N; ++i) {
        for (size_t j = 0; j < KYBER_K; ++j) {
            for (size_t k = 0; k < 2 * KYBER_ETA1; ++k) {
                V[i].s_ddeg_sh[j][k] = s_reduced_ddeg_shares[j][k].share_y[i];
                V[i].e_ddeg_sh[j][k] = e_reduced_ddeg_shares[j][k].share_y[i];
                V[i].s_zero_sh[j][k] = s_zero_2ddeg_shares[j][k].share_y[i];
                V[i].e_zero_sh[j][k] = e_zero_2ddeg_shares[j][k].share_y[i];
            }
        }
    }

    // build challenge:
    uint8_t ch_seeds[MPCITH_N][KYBER_SYMBYTES];
    // extract the views from challenge:
    for (size_t i = 0; i < MPCITH_N; ++i) {
        uint8_t view[KYBER_SYMBYTES + (((6 + 4 * KYBER_ETA1 * 2) * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1))) * sizeof(uint16_t)] = {0};
        // copy i-th view to view[]
        memcpy(view,
               V[i].comm,
               KYBER_SYMBYTES);
        memcpy(view + KYBER_SYMBYTES,
               V[i].s_sh,
               KYBER_K * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + KYBER_K * sizeof(uint16_t),
               V[i].e_sh,
               KYBER_K * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + (2 * KYBER_K) * sizeof(uint16_t),
               V[i].f_sh,
               (MPCITH_K + MPCITH_V + 1) * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + (2 * KYBER_K + (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t),
               V[i].Tf_sh,
               (MPCITH_K + MPCITH_V + 1) * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + (2 * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t),
               V[i].beta,
               KYBER_K * sizeof(uint16_t));
        memcpy(view+ KYBER_SYMBYTES + (3 * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t),
               V[i].gamma,
               KYBER_K * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + (4 * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t),
               V[i].sr_sh,
               KYBER_K * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + (5 * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t),
               V[i].er_sh,
               KYBER_K * sizeof(uint16_t));

        for (size_t j = 0; j < KYBER_K; ++j) {
            size_t curr_offset = KYBER_SYMBYTES + (6 * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1) + 4 * (j * 2 * KYBER_ETA1)) * sizeof(uint16_t);
            memcpy(view + curr_offset,
                   V[i].s_ddeg_sh[j],
                   2 * KYBER_ETA1 * sizeof(uint16_t));
            memcpy(view + curr_offset + 1 * (2 * KYBER_ETA1) * sizeof(uint16_t),
                   V[i].e_ddeg_sh[j],
                   2 * KYBER_ETA1 * sizeof(uint16_t));
            memcpy(view + curr_offset + 2 * (2 * KYBER_ETA1) * sizeof(uint16_t),
                   V[i].s_zero_sh[j],
                   2 * KYBER_ETA1 * sizeof(uint16_t));
            memcpy(view + curr_offset + 3 * (2 * KYBER_ETA1) * sizeof(uint16_t),
                   V[i].e_zero_sh[j],
                   2 * KYBER_ETA1 * sizeof(uint16_t));
        }
        sha3_256(ch_seeds[i], view, sizeof(view));
    }
    uint8_t merged_ch[MPCITH_N * KYBER_SYMBYTES], ch[KYBER_SYMBYTES];
    for (size_t i = 0; i < MPCITH_N; ++i) {
        memcpy(merged_ch + i * KYBER_SYMBYTES, ch_seeds[i], KYBER_SYMBYTES);
    }
    sha3_256(ch, merged_ch, MPCITH_N * KYBER_SYMBYTES);

    // selected open:
    uint8_t I_[2 * MPCITH_T];
    kyber_shake256_prf(I_, 2 * MPCITH_T, ch, 1);
    // merge I_ into I and mod 1454
    for (size_t i = 0; i < MPCITH_T; ++i) {
        pi->I[i] = ((I_[2 * i] << 8) | I_[2 * i + 1]) % MPCITH_N;
    }
    // remove duplicate
    for (size_t i = 1; i < MPCITH_T; i++) {
        size_t j = 0;
        uint16_t inc = 0;
        bool is_dup;
        do {
            is_dup = false;
            for (j = 0; j < i; j++) {
                if ((pi->I[i] + inc) % MPCITH_N == pi->I[j]) {
                    is_dup = true;
                    inc = inc + 1;
                    break;
                }
            }
        } while (is_dup);
        pi->I[i] = (pi->I[i] + inc) % MPCITH_N;
    }

    /**
     *  build proof:
     */
    // compute rest index of parties
    bool in_I[MPCITH_N] = {false};
    uint16_t rest_I[MPCITH_N - MPCITH_T];
    for (size_t i = 0; i < MPCITH_T; ++i) {
        in_I[pi->I[i]] = true;
    }
    for (size_t i = 0, j = 0; i < MPCITH_N; ++i) {
        if (!in_I[i]) {
            rest_I[j] = i;
            ++j;
        }
    }

    // views in set I
    for (size_t i = 0; i < MPCITH_T; ++i) {
        for (size_t j = 0; j < MPCITH_K + MPCITH_V + 1; ++j) {
            pi->f_shares[i][j] = V[pi->I[i]].f_sh[j];
            pi->NTT_f_shares[i][j] = V[pi->I[i]].Tf_sh[j];
        }
        for (size_t j = 0; j < KYBER_K; ++j) {
            pi->s_shares[i][j] = V[pi->I[i]].s_sh[j];
            pi->e_shares[i][j] = V[pi->I[i]].e_sh[j];
            pi->NTT_s_shares[i][j] = ntt_s_shares[j].share_y[pi->I[i]];
            pi->NTT_e_shares[i][j] = ntt_e_shares[j].share_y[pi->I[i]];
            pi->NTT_Ar_shares[i][j] = ntt_Ar_shares[j].share_y[pi->I[i]];
            pi->NTT_As_shares[i][j] = ntt_As_shares[j].share_y[pi->I[i]];
            for (size_t k = 0; k < 2 * KYBER_ETA1 + 1; ++k) {
                pi->s_sub_eta_shares[i][j][k] = s_sub_eta_shares[j][k].share_y[pi->I[i]];
                pi->e_sub_eta_shares[i][j][k] = e_sub_eta_shares[j][k].share_y[pi->I[i]];
            }
            for (size_t k = 0; k < 2 * KYBER_ETA1; ++k) {
                pi->z_s_ddeg_shares[i][j][k] = V[pi->I[i]].s_ddeg_sh[j][k];
                pi->z_e_ddeg_shares[i][j][k] = V[pi->I[i]].e_ddeg_sh[j][k];
            }
        }
    }
    // views in set [N - T]
    for (size_t i = 0; i < MPCITH_N - MPCITH_T; ++i) {
        for (size_t j = 0; j < MPCITH_K; ++j) {
            pi->beta_shares[i][j] = V[rest_I[i]].beta[j];
            pi->gamma_shares[i][j] = V[rest_I[i]].gamma[j];
        }
        for (size_t j = 0; j < KYBER_K; ++j) {
            pi->sr_shares[i][j] = V[rest_I[i]].sr_sh[j];
            pi->er_shares[i][j] = V[rest_I[i]].er_sh[j];
            pi->t_shares[i][j] = ntt_t_shares[j].share_y[rest_I[i]];
            for (size_t k = 0; k < 2 * KYBER_ETA1 + 1; ++k) {
                pi->s_eta_shares[i][j][k] = eta_share->s_eta_shares[j][k].share_y[rest_I[i]];
                pi->e_eta_shares[i][j][k] = eta_share->e_eta_shares[j][k].share_y[rest_I[i]];
            }
            for (size_t k = 0; k < 2 * KYBER_ETA1; ++k) {
                pi->u_s_2ddeg_shares[i][j][k] = V[rest_I[i]].s_zero_sh[j][k];
                pi->u_e_2ddeg_shares[i][j][k] = V[rest_I[i]].e_zero_sh[j][k];
            }

        }
        memcpy(pi->Tcomm[i], V[rest_I[i]].comm, KYBER_SYMBYTES);
        memcpy(pi->comm[i], ch_seeds[rest_I[i]], KYBER_SYMBYTES);
    }
}

void encode_mpcith_proof(uint8_t *buf,
                         const mpcith_proof *pi) {
    memcpy(buf, pi, sizeof(mpcith_proof));
}

void decode_mpcith_proof(mpcith_proof *pi, const uint8_t *buf) {
    memcpy(pi->f_shares,
           buf,
           (MPCITH_T * (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t));
    memcpy(pi->NTT_f_shares,
           buf + (MPCITH_T * (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t),
           (MPCITH_T * (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t));
    size_t offset = MPCITH_T * (MPCITH_K + MPCITH_V + 1) * 2 * sizeof(uint16_t);
    memcpy(pi->beta_shares,
           buf + offset,
           (MPCITH_N - MPCITH_T) * MPCITH_K * sizeof(uint16_t));
    memcpy(pi->gamma_shares,
           buf + offset + (MPCITH_N - MPCITH_T) * MPCITH_K * sizeof(uint16_t),
           (MPCITH_N - MPCITH_T) * MPCITH_K * sizeof(uint16_t));
    offset += (MPCITH_N - MPCITH_T) * MPCITH_K * 2 * sizeof(uint16_t);
    memcpy(pi->Tcomm,
           buf + offset,
           (MPCITH_N - MPCITH_T) * KYBER_SYMBYTES * sizeof(uint8_t));
    offset += (MPCITH_N - MPCITH_T) * KYBER_SYMBYTES * sizeof(uint8_t);
    memcpy(pi->I,
           buf + offset,
           MPCITH_T * sizeof(uint16_t));
    offset += MPCITH_T * sizeof(uint16_t);
    memcpy(pi->s_shares,
           buf + offset,
           MPCITH_T * KYBER_K * sizeof(uint16_t));
    memcpy(pi->e_shares,
           buf + offset + MPCITH_T * KYBER_K * sizeof(uint16_t),
           MPCITH_T * KYBER_K * sizeof(uint16_t));
    offset += MPCITH_T * KYBER_K * 2 * sizeof(uint16_t);
    memcpy(pi->t_shares,
           buf + offset,
           (MPCITH_N - MPCITH_T) * KYBER_K * sizeof(uint16_t));
    offset += (MPCITH_N - MPCITH_T) * KYBER_K * sizeof(uint16_t);
    memcpy(pi->NTT_s_shares,
           buf + offset,
           MPCITH_T * KYBER_K * sizeof(uint16_t));
    memcpy(pi->NTT_e_shares,
           buf + offset + MPCITH_T * KYBER_K * sizeof(uint16_t),
           MPCITH_T * KYBER_K * sizeof(uint16_t));
    memcpy(pi->NTT_Ar_shares,
           buf + offset + MPCITH_T * KYBER_K * 2 * sizeof(uint16_t),
           MPCITH_T * KYBER_K * sizeof(uint16_t));
    memcpy(pi->NTT_As_shares,
           buf + offset + MPCITH_T * KYBER_K * 3 * sizeof(uint16_t),
           MPCITH_T * KYBER_K * sizeof(uint16_t));
    offset += MPCITH_T * KYBER_K * 4 * sizeof(uint16_t);
    memcpy(pi->sr_shares,
           buf + offset,
           (MPCITH_N - MPCITH_T) * KYBER_K * sizeof(uint16_t));
    memcpy(pi->er_shares,
           buf + offset + (MPCITH_N - MPCITH_T) * KYBER_K * sizeof(uint16_t),
           (MPCITH_N - MPCITH_T) * KYBER_K * sizeof(uint16_t));
    offset += (MPCITH_N - MPCITH_T) * KYBER_K * 2 * sizeof(uint16_t);
    memcpy(pi->s_eta_shares,
           buf + offset,
           (MPCITH_N - MPCITH_T) * KYBER_K * (KYBER_ETA1 * 2 + 1) * sizeof(uint16_t));
    memcpy(pi->e_eta_shares,
           buf + offset + (MPCITH_N - MPCITH_T) * KYBER_K * (KYBER_ETA1 * 2 + 1) * sizeof(uint16_t),
           (MPCITH_N - MPCITH_T) * KYBER_K * (KYBER_ETA1 * 2 + 1) * sizeof(uint16_t));
    offset += (MPCITH_N - MPCITH_T) * KYBER_K * (KYBER_ETA1 * 2 + 1) * 2 * sizeof(uint16_t);
    memcpy(pi->s_sub_eta_shares,
           buf + offset,
           MPCITH_T * KYBER_K * (KYBER_ETA1 * 2 + 1) * sizeof(uint16_t));
    memcpy(pi->e_sub_eta_shares,
           buf + offset + MPCITH_T * KYBER_K * (KYBER_ETA1 * 2 + 1) * sizeof(uint16_t),
           MPCITH_T * KYBER_K * (KYBER_ETA1 * 2 + 1) * sizeof(uint16_t));
    offset += MPCITH_T * KYBER_K * (KYBER_ETA1 * 2 + 1) * 2 * sizeof(uint16_t);
    memcpy(pi->z_s_ddeg_shares,
           buf + offset,
           MPCITH_T * KYBER_K * (KYBER_ETA1 * 2) * sizeof(uint16_t));
    memcpy(pi->z_e_ddeg_shares,
           buf + offset + MPCITH_T * KYBER_K * (KYBER_ETA1 * 2) * sizeof(uint16_t),
           MPCITH_T * KYBER_K * (KYBER_ETA1 * 2) * sizeof(uint16_t));
    offset += MPCITH_T * KYBER_K * (KYBER_ETA1 * 2) * 2 * sizeof(uint16_t);
    memcpy(pi->u_s_2ddeg_shares,
           buf + offset,
           (MPCITH_N - MPCITH_T) * KYBER_K * (KYBER_ETA1 * 2) * sizeof(uint16_t));
    memcpy(pi->u_e_2ddeg_shares,
           buf + offset + (MPCITH_N - MPCITH_T) * KYBER_K * (KYBER_ETA1 * 2) * sizeof(uint16_t),
           (MPCITH_N - MPCITH_T) * KYBER_K * (KYBER_ETA1 * 2) * sizeof(uint16_t));
    offset += (MPCITH_N - MPCITH_T) * KYBER_K * (KYBER_ETA1 * 2) * 2 * sizeof(uint16_t);
    memcpy(pi->comm,
           buf + offset,
           (MPCITH_N - MPCITH_T) * KYBER_SYMBYTES * sizeof(uint8_t));
}