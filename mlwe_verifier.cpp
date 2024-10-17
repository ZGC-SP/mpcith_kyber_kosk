
#include "mlwe_verifier.hpp"

bool verify(const mpcith_proof *pi, const mlwe_inst *mlwe) {
    // prepare the NTL module
    ZZ_p::init(ZZ(KYBER_Q));

    // compute the set I and rest set:
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

    // recompute: Tcomm[i] = Hash(s_sh_i, e_sh_i, f_sh_i, T_sh_i)
    uint8_t Tcomm_rec[MPCITH_N][KYBER_SYMBYTES];
    for (size_t i = 0; i < MPCITH_T; i++) {
        uint16_t ctxt[(KYBER_K + MPCITH_K + MPCITH_V + 1) * 2];
        for (size_t j = 0; j < KYBER_K; ++j) {
            ctxt[j] = pi->s_shares[i][j];
            ctxt[j + KYBER_K] = pi->e_shares[i][j];

        }
        for (size_t j = 0; j < MPCITH_K + MPCITH_V + 1; ++j) {
            ctxt[j + (KYBER_K * 2)] = pi->f_shares[i][j];
            ctxt[j + (KYBER_K * 2) + (MPCITH_K + MPCITH_V + 1)] = pi->NTT_f_shares[i][j];
        }
        sha3_256(Tcomm_rec[pi->I[i]], (const uint8_t *) ctxt, sizeof(ctxt));
    }
    for (size_t i = 0; i < MPCITH_N - MPCITH_T; i++) {
        memcpy(Tcomm_rec[rest_I[i]], pi->Tcomm[i], KYBER_SYMBYTES);
    }

    // recompute alpha:
    uint8_t alpha_[sizeof(uint16_t) * (MPCITH_K + MPCITH_V)];
    uint8_t Tcomm[MPCITH_N * KYBER_SYMBYTES], digest[KYBER_SYMBYTES];
    for (size_t i = 0; i < MPCITH_N; ++i) {
        memcpy(Tcomm + i * KYBER_SYMBYTES, Tcomm_rec[i], KYBER_SYMBYTES);
    }
    sha3_256(digest, Tcomm, MPCITH_N * KYBER_SYMBYTES);
    kyber_shake256_prf(alpha_, sizeof(uint16_t) * (MPCITH_K + MPCITH_V), digest, 1);

    uint16_t alpha[MPCITH_K + MPCITH_V];
    // encode uint8 to GF(3329)
    for (size_t i = 0; i < MPCITH_K + MPCITH_V; ++i) {
        alpha[i] = ((alpha_[2 * i] << 8) | alpha_[2 * i + 1]) % KYBER_Q;
    }

    // recompute (MPCITH_T) beta, gamma via poly eval
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

    uint16_t beta[MPCITH_N][MPCITH_K], gamma[MPCITH_N][MPCITH_K];
    for (size_t i = 0; i < MPCITH_T; ++i) {
        for (size_t j = 0; j < MPCITH_K; ++j) {
            beta[pi->I[i]][j] = 0;
            for (size_t k = 0; k < MPCITH_K + MPCITH_V + 1; ++k) {
                if (k == 0) {
                    beta[pi->I[i]][j] = pi->f_shares[i][k];
                } else {
                    beta[pi->I[i]][j] = gf3329_add(beta[pi->I[i]][j],
                                                   gf3329_mul(pow_table[j][k], pi->f_shares[i][k]));
                }
            }
            gamma[pi->I[i]][j] = 0;
            for (size_t k = 0; k < MPCITH_K + MPCITH_V + 1; ++k) {
                if (k == 0) {
                    gamma[pi->I[i]][j] = pi->NTT_f_shares[i][k];
                } else {
                    gamma[pi->I[i]][j] = gf3329_add(gamma[pi->I[i]][j],
                                                    gf3329_mul(pow_table[j][k], pi->NTT_f_shares[i][k]));
                }
            }
        }
    }
    // rest beta and gamma copy from pi
    for (size_t i = 0; i < MPCITH_N - MPCITH_T; ++i) {
        for (size_t j = 0; j < MPCITH_K; ++j) {
            beta[rest_I[i]][j] = pi->beta_shares[i][j];
            gamma[rest_I[i]][j] = pi->gamma_shares[i][j];
        }
    }
    share_vec beta_shares[MPCITH_K], gamma_shares[MPCITH_K];
    secret_vec beta_secret[MPCITH_K], gamma_secret[MPCITH_K];
    for (size_t i = 0; i < MPCITH_K; ++i) {
        for (size_t j = 0; j < MPCITH_N; ++j) {
            beta_shares[i].share_x[j] = j + KYBER_N;
            beta_shares[i].share_y[j] = beta[j][i];
            gamma_shares[i].share_x[j] = j + KYBER_N;
            gamma_shares[i].share_y[j] = gamma[j][i];
        }
        recon_secrets_ddeg(&beta_secret[i], &beta_shares[i]);
        recon_secrets_ddeg(&gamma_secret[i], &gamma_shares[i]);
    }

    // check beta and gamma relation: gamma = NTT(beta)
    for (size_t i = 0; i < MPCITH_K; ++i) {
        poly beta_poly, gamma_poly;
        for (size_t j = 0; j < KYBER_N; ++j) {
            beta_poly.coeffs[j] = decode_from_gf3329(beta_secret[i].secret[j]);
            gamma_poly.coeffs[j] = decode_from_gf3329(gamma_secret[i].secret[j]);
        }
        poly_ntt(&beta_poly);
        for (size_t j = 0; j < KYBER_N; ++j) {
            if (beta_poly.coeffs[j] != gamma_poly.coeffs[j]) {
                printf("Check failed for beta[%zu] and gamma[%zu] relation.\n", i, i);
                return false;
            }
        }
    }
    // check polynomial constraints
    for (size_t i = 0; i < MPCITH_K; ++i) {
        for (size_t j = 0, k = 0; j < MPCITH_N; ++j) {
            if (!in_I[j]) {
                if (beta[j][i] != pi->beta_shares[k][i]) {
                    printf("Check failed for beta share: (rec)%d, (org)%d at %zu.\n",
                           beta[j][i], pi->beta_shares[k][i], j + KYBER_N);
                    return false;
                }
                if (gamma[j][i] != pi->gamma_shares[k][i]) {
                    printf("Check failed for gamma share: (rec)%d, (org)%d at %zu.\n",
                           gamma[j][i], pi->gamma_shares[k][i], j + KYBER_N);
                    return false;
                }
                k += 1;
            }
        }
    }

    /**
     * Verify the A * s + e relation.
     */
    // recompute r and NTT_r
    uint16_t r[MPCITH_T][MPCITH_V], NTT_r[MPCITH_T][MPCITH_V];
    for (size_t i = 0; i < MPCITH_T; ++i) {
        for (size_t j = MPCITH_K; j < MPCITH_V + MPCITH_K; ++j) {
            r[i][j - MPCITH_K] = 0;
            for (size_t k = 0; k < MPCITH_K + MPCITH_V + 1; ++k) {
                if (k == 0) {
                    r[i][j - MPCITH_K] = pi->f_shares[i][k + MPCITH_K + 1];
                } else {
                    r[i][j - MPCITH_K] = gf3329_add(r[i][j - MPCITH_K],
                                                    gf3329_mul(pow_table[j][k], pi->f_shares[i][k]));
                }
            }
            NTT_r[i][j - MPCITH_K] = 0;
            for (size_t k = 0; k < MPCITH_K + MPCITH_V + 1; ++k) {
                if (k == 0) {
                    NTT_r[i][j - MPCITH_K] = pi->NTT_f_shares[i][k + MPCITH_K + 1];
                } else {
                    NTT_r[i][j - MPCITH_K] = gf3329_add(NTT_r[i][j - MPCITH_K],
                                                        gf3329_mul(pow_table[j][k], pi->NTT_f_shares[i][k]));
                }
            }
        }
    }

    // recon s + r1 and e + r2
    uint16_t sr_yval[KYBER_K][DEG_D + 1], er_yval[KYBER_K][DEG_D + 1];
    uint16_t sr_rnd[KYBER_K][DEG_D + 1], er_rnd[KYBER_K][DEG_D + 1];
    share_vec s_add_r_shares[KYBER_K], e_add_r_shares[KYBER_K];
    polyvec s_add_r, e_add_r;
    for (size_t i = 0; i < KYBER_K; ++i) {
        uint16_t sr_x[DEG_D + 1], sr_y[DEG_D + 1];
        uint16_t er_x[DEG_D + 1], er_y[DEG_D + 1];
        for (size_t j = 0; j < DEG_D + 1; ++j) {
            sr_x[j] = rest_I[j] + KYBER_N;
            sr_y[j] = pi->sr_shares[j][i];

            er_x[j] = rest_I[j] + KYBER_N;
            er_y[j] = pi->er_shares[j][i];
        }

        Vec<ZZ_p> ntl_sr_x, ntl_sr_y, ntl_er_x, ntl_er_y;
        ntl_sr_x.SetLength(DEG_D + 1);
        ntl_sr_y.SetLength(DEG_D + 1);
        ntl_er_x.SetLength(DEG_D + 1);
        ntl_er_y.SetLength(DEG_D + 1);
        for (size_t j = 0; j < DEG_D + 1; ++j) {
            ntl_sr_x[j] = sr_x[j];
            ntl_sr_y[j] = sr_y[j];
            ntl_er_x[j] = er_x[j];
            ntl_er_y[j] = er_y[j];
        }

        ZZ_pX sr_poly, er_poly;
        interpolate(sr_poly, ntl_sr_x, ntl_sr_y);
        interpolate(er_poly, ntl_er_x, ntl_er_y);
        // eval the first KYBER_N points for s + r and e + r:
        Vec<ZZ_p> sr_y_ntl, er_y_ntl;
        sr_y_ntl.SetLength(KYBER_N);
        er_y_ntl.SetLength(KYBER_N);
        for (size_t j = 0; j < KYBER_N; ++j) {
            sr_y_ntl[j] = eval(sr_poly, ZZ_p(j));
            er_y_ntl[j] = eval(er_poly, ZZ_p(j));
            // encode to polyvec
            s_add_r.vec[i].coeffs[j] = gf3329_int_to_int16(conv<int>(sr_y_ntl[j]));
            e_add_r.vec[i].coeffs[j] = gf3329_int_to_int16(conv<int>(er_y_ntl[j]));
            sr_yval[i][j] = gf3329_from_int(conv<int>(sr_y_ntl[j]));
            er_yval[i][j] = gf3329_from_int(conv<int>(er_y_ntl[j]));
        }
        // rest points to hold degree
        for (size_t j = KYBER_N; j < DEG_D + 1; ++j) {
            sr_yval[i][j] = gf3329_from_int(conv<int>(eval(sr_poly, ZZ_p(j))));
            er_yval[i][j] = gf3329_from_int(conv<int>(eval(er_poly, ZZ_p(j))));
        }

        // reshare the polynomial shares

        recompute_share_secrets_ddeg(&s_add_r_shares[i], (const uint16_t *) sr_yval[i]);
        recompute_share_secrets_ddeg(&e_add_r_shares[i], (const uint16_t *) er_yval[i]);
        // extract randomness
        for (size_t j = KYBER_N; j < DEG_D + 1; ++j) {
            sr_rnd[i][j] = s_add_r_shares[i].share_y[j - KYBER_N];
            er_rnd[i][j] = e_add_r_shares[i].share_y[j - KYBER_N];
        }
        // validate the polynomial constraints of s + r and e + r.
        for (size_t j = 0, p_i = 0; j < MPCITH_N; ++j) {
            if (!in_I[j]) {
                if (s_add_r_shares[i].share_y[j] != pi->sr_shares[p_i][i]) {
                    printf("s + r share error: (rec)%d, (org)%d at %zu.\n",
                           s_add_r_shares[i].share_y[j], pi->sr_shares[p_i][i], j + KYBER_N);
                    return false;
                }
                if (e_add_r_shares[i].share_y[j] != pi->er_shares[p_i][i]) {
                    printf("e + r share error: (rec)%d, (org)%d at %zu.\n",
                           e_add_r_shares[i].share_y[j], pi->er_shares[p_i][i], j + KYBER_N);
                    return false;
                }
                p_i += 1;
            }
        }
    }

    uint16_t sr_shares_view[MPCITH_T][KYBER_K], er_shares_view[MPCITH_T][KYBER_K];
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < MPCITH_T; ++j) {
            sr_shares_view[j][i] = s_add_r_shares[i].share_y[pi->I[j]];
            er_shares_view[j][i] = e_add_r_shares[i].share_y[pi->I[j]];
        }
    }

    polyvec_ntt(&s_add_r);
    polyvec_ntt(&e_add_r);
    // set rnd
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_N; ++j) {
            sr_rnd[i][j] = encode_to_gf3329(s_add_r.vec[i].coeffs[j]);
            er_rnd[i][j] = encode_to_gf3329(e_add_r.vec[i].coeffs[j]);
        }
    }

    share_vec ntt_sr_shares[KYBER_K], ntt_er_shares[KYBER_K];
    for (size_t i = 0; i < KYBER_K; ++i) {
        recompute_share_secrets_ddeg(&ntt_sr_shares[i], sr_rnd[i]);
        recompute_share_secrets_ddeg(&ntt_er_shares[i], er_rnd[i]);
    }
    // validate opened shares of NTT(s) = NTT(s + r) - NTT(r) and also for e
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < MPCITH_T; ++j) {
            if (pi->NTT_s_shares[j][i] != gf3329_sub(ntt_sr_shares[i].share_y[pi->I[j]], NTT_r[j][i])) {
                printf("Check failed for NTT(s[%zu]) at view %d.\n", i, pi->I[j]);
                return false;
            }
            if (pi->NTT_e_shares[j][i] != gf3329_sub(ntt_er_shares[i].share_y[pi->I[j]], NTT_r[j][i + KYBER_K])) {
                printf("Check failed for NTT(e[%zu]) at view %d.\n", i, pi->I[j]);
                return false;
            }
        }
    }


    polyvec ntt_Asr; share_vec ntt_Asr_shares[KYBER_K];
    for (size_t i = 0; i < KYBER_K; ++i) {
        polyvec_basemul_acc_montgomery(&ntt_Asr.vec[i], &mlwe->A[i], &s_add_r);
        poly_tomont(&ntt_Asr.vec[i]);
    }
    for (size_t i = 0; i < KYBER_K; ++i) {
        uint16_t ntt_Asr_rnd[DEG_D + 1];
        for (size_t j = 0; j < KYBER_N; ++j) {
            ntt_Asr_rnd[j] = encode_to_gf3329(ntt_Asr.vec[i].coeffs[j]);
        }
        for (size_t j = KYBER_N; j < DEG_D + 1; ++j) {
            ntt_Asr_rnd[j] = sr_rnd[i][j];
        }
        recompute_share_secrets_ddeg(&ntt_Asr_shares[i], ntt_Asr_rnd);
    }

    // Check: NTT(A*(s+r)) = NTT(A*s) + NTT(A*r)
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < MPCITH_T; ++j) {
            if (ntt_Asr_shares[i].share_y[pi->I[j]] != gf3329_add(pi->NTT_As_shares[j][i], pi->NTT_Ar_shares[j][i])) {
                printf("Check failed for NTT(A*(s[%zu]+r)) at view %d.\n",
                       i, pi->I[j]);
                return false;
            }
        }
    }

    // check A * s + e = t:
    // recounstruct the t:
    secret_vec t_secret[KYBER_K]; share_vec t_shares[KYBER_K];
    polyvec t; uint16_t t_yval[KYBER_K][DEG_D + 1];

    for (size_t i = 0; i < KYBER_K; ++i) {
        uint16_t t_x[DEG_D + 1], t_y[DEG_D + 1];
        for (size_t j = 0; j < DEG_D + 1; ++j) {
            t_x[j] = rest_I[j] + KYBER_N;
            t_y[j] = pi->t_shares[j][i];
        }

        Vec<ZZ_p> ntl_t_x, ntl_t_y;
        ntl_t_x.SetLength(DEG_D + 1);
        ntl_t_y.SetLength(DEG_D + 1);
        for (size_t j = 0; j < DEG_D + 1; ++j) {
            ntl_t_x[j] = t_x[j];
            ntl_t_y[j] = t_y[j];
        }

        clock_t start_interpolate_ntl = clock();

        ZZ_pX t_poly;
        interpolate(t_poly, ntl_t_x, ntl_t_y);
        // eval the first KYBER_N points for t:
        Vec<ZZ_p> t_y_ntl;
        t_y_ntl.SetLength(KYBER_N);
        for (size_t j = 0; j < KYBER_N; ++j) {
            t_y_ntl[j] = eval(t_poly, ZZ_p(j));
            // encode to polyvec
            t.vec[i].coeffs[j] = gf3329_int_to_int16(conv<int>(t_y_ntl[j]));
            t_yval[i][j] = gf3329_from_int(conv<int>(t_y_ntl[j]));
        }
        // eval the rest points for t:
        for (size_t j = KYBER_N; j < DEG_D + 1; ++j) {
            t_yval[i][j] = gf3329_from_int(conv<int>(eval(t_poly, ZZ_p(j))));
        }
        recompute_share_secrets_ddeg(&t_shares[i], (const uint16_t *) t_yval[i]);
    }

    polyvec_reduce(&t);
    // validate the recomputed t is equal to t from pk or not
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_N; ++j) {
            if (encode_to_gf3329(t.vec[i].coeffs[j]) != encode_to_gf3329(mlwe->t.vec[i].coeffs[j])) {
                printf("Check failed for t[%zu] at %zu.\n", i, j);
                return false;
            }
        }
    }
    // validate the relation of A * NTT(s) + NTT(e) = NTT(t) holds or not
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < MPCITH_T; ++j) {
            uint16_t As_ = pi->NTT_As_shares[j][i];
            uint16_t e_ = pi->NTT_e_shares[j][i];
            uint16_t t_ = t_shares[i].share_y[pi->I[j]];
            if (t_ != gf3329_add(As_, e_)) {
                printf("Check failed for t = A*s + e at view %d.\n",
                       i, pi->I[j]);
                return false;
            }
        }
    }

    /**
     * Verify the s and e range over [-eta, eta].
     */
    // recon the eta range shares: [-eta , eta]
    uint16_t s_eta_yval[KYBER_K][KYBER_ETA1 * 2 + 1][DEG_D + 1];
    uint16_t e_eta_yval[KYBER_K][KYBER_ETA1 * 2 + 1][DEG_D + 1];
    share_vec s_eta_shares[KYBER_K][KYBER_ETA1 * 2 + 1];
    share_vec e_eta_shares[KYBER_K][KYBER_ETA1 * 2 + 1];
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_ETA1 * 2 + 1; ++j) {
            uint16_t s_eta_x[DEG_D + 1], s_eta_y[DEG_D + 1];
            uint16_t e_eta_x[DEG_D + 1], e_eta_y[DEG_D + 1];
            for (size_t k = 0; k < DEG_D + 1; ++k) {
                s_eta_x[k] = rest_I[k] + KYBER_N;
                e_eta_x[k] = rest_I[k] + KYBER_N;
                s_eta_y[k] = pi->s_eta_shares[k][i][j];
                e_eta_y[k] = pi->e_eta_shares[k][i][j];
            }

            Vec<ZZ_p> ntl_s_eta_x, ntl_s_eta_y, ntl_e_eta_x, ntl_e_eta_y;
            ntl_s_eta_x.SetLength(DEG_D + 1);
            ntl_s_eta_y.SetLength(DEG_D + 1);
            ntl_e_eta_x.SetLength(DEG_D + 1);
            ntl_e_eta_y.SetLength(DEG_D + 1);
            for (size_t k = 0; k < DEG_D + 1; ++k) {
                ntl_s_eta_x[k] = s_eta_x[k];
                ntl_e_eta_x[k] = e_eta_x[k];

                ntl_s_eta_y[k] = s_eta_y[k];
                ntl_e_eta_y[k] = e_eta_y[k];
            }
            ZZ_pX s_eta_poly, e_eta_poly;
            interpolate(s_eta_poly, ntl_s_eta_x, ntl_s_eta_y);
            interpolate(e_eta_poly, ntl_e_eta_x, ntl_e_eta_y);
            // eval the first KYBER_N points for eta shares for s / e:
            Vec<ZZ_p> s_eta_y_ntl, e_eta_y_ntl;
            e_eta_y_ntl.SetLength(KYBER_N), s_eta_y_ntl.SetLength(KYBER_N);
            for (size_t k = 0; k < KYBER_N; ++k) {
                s_eta_y_ntl[k] = eval(s_eta_poly, ZZ_p(k));
                e_eta_y_ntl[k] = eval(e_eta_poly, ZZ_p(k));
                uint16_t current_eta = gf3329_sub(j, KYBER_ETA1);
                // validate the eta shares for s and e
                if (current_eta != gf3329_from_int(conv<int>(s_eta_y_ntl[k]))) {
                    printf("Check failed for s_eta[%zu] at %zu (%d != %d).\n",
                           i, k, current_eta, conv<int>(s_eta_y_ntl[k]));
                    return false;
                }
                if (current_eta != gf3329_from_int(conv<int>(e_eta_y_ntl[k]))) {
                    printf("Check failed for e_eta[%zu] at %zu (%d != %d).\n",
                           i, k, current_eta, conv<int>(e_eta_y_ntl[k]));
                    return false;
                }
            }
            // eval the rest points for t:
            for (size_t k = 0; k < KYBER_N; ++k) {
                uint16_t current_eta = gf3329_sub(j, KYBER_ETA1);
                s_eta_yval[i][j][k] = current_eta;
                e_eta_yval[i][j][k] = current_eta;
            }
            for (size_t k = KYBER_N; k < DEG_D + 1; ++k) {
                s_eta_yval[i][j][k] = gf3329_from_int(conv<int>(eval(s_eta_poly, ZZ_p(k))));
                e_eta_yval[i][j][k] = gf3329_from_int(conv<int>(eval(e_eta_poly, ZZ_p(k))));
            }
            recompute_share_secrets_ddeg(&s_eta_shares[i][j], (const uint16_t *) s_eta_yval[i][j]);
            recompute_share_secrets_ddeg(&e_eta_shares[i][j], (const uint16_t *) e_eta_yval[i][j]);
        }
    }

    // recompute (s - eta[i]) and (e - eta[i]) for opened views
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < MPCITH_T; ++j) {
            for (size_t k = 0; k < KYBER_ETA1 * 2 + 1; ++k) {
                // check [s] - [eta] =? [s - eta] from recomputed shares.
                if (pi->s_sub_eta_shares[j][i][k] != gf3329_sub(pi->s_shares[j][i],
                                                                s_eta_shares[i][k].share_y[pi->I[j]])) {
                    printf("Check failed for s - eta at view %d, eta[%d][%d][%d].\n",
                           pi->I[j],
                           i,j,k);
                    return false;
                }
                // check [e] - [eta] =? [e - eta] from recomputed shares.
                if (pi->e_sub_eta_shares[j][i][k] != gf3329_sub(pi->e_shares[j][i],
                                                                e_eta_shares[i][k].share_y[pi->I[j]])) {
                    printf("Check failed for e - eta at view %d, eta %zu.\n", pi->I[j], j - KYBER_ETA1);
                    return false;
                }
            }
        }
    }

    // eval the [z_2d] = [x] * [y] and compute [u_2d] = [z_2d] - z[d]
    share_vec s_u_shares[KYBER_K][KYBER_ETA1 * 2], e_u_shares[KYBER_K][KYBER_ETA1 * 2];
    uint16_t s_z2d_sh[KYBER_K][2 * KYBER_ETA1][MPCITH_T], e_z2d_sh[KYBER_K][2 * KYBER_ETA1][MPCITH_T];
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < 2 * KYBER_ETA1; ++j) {
            for (size_t k = 0; k < MPCITH_T; ++k) {
                if (j == 0) {
                    s_z2d_sh[i][j][k] = gf3329_mul(pi->s_sub_eta_shares[k][i][j],
                                                   pi->s_sub_eta_shares[k][i][j + 1]);
                    e_z2d_sh[i][j][k] = gf3329_mul(pi->e_sub_eta_shares[k][i][j],
                                                   pi->e_sub_eta_shares[k][i][j + 1]);
                } else {
                    s_z2d_sh[i][j][k] = gf3329_mul(pi->z_s_ddeg_shares[k][i][j - 1],
                                                   pi->s_sub_eta_shares[k][i][j + 1]);
                    e_z2d_sh[i][j][k] = gf3329_mul(pi->z_e_ddeg_shares[k][i][j - 1],
                                                   pi->e_sub_eta_shares[k][i][j + 1]);
                }
                // set shares for u_2d
                s_u_shares[i][j].share_x[pi->I[k]] = pi->I[k] + KYBER_N;
                s_u_shares[i][j].share_y[pi->I[k]] = gf3329_sub(s_z2d_sh[i][j][k],
                                                                pi->z_s_ddeg_shares[k][i][j]);

                e_u_shares[i][j].share_x[pi->I[k]] = pi->I[k] + KYBER_N;
                e_u_shares[i][j].share_y[pi->I[k]] = gf3329_sub(e_z2d_sh[i][j][k],
                                                                pi->z_e_ddeg_shares[k][i][j]);
            }
        }
    }

    // recon the [u] = [z_2d] - z[d] = [0]
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < KYBER_ETA1 * 2; ++j) {
            // interpolate the u2d from N
            uint16_t s_u_x[DEG_2D + 1], s_u_y[DEG_2D + 1];
            uint16_t e_u_x[DEG_2D + 1], e_u_y[DEG_2D + 1];
            for (size_t k = 0; k < DEG_2D + 1; ++k) {
                s_u_x[k] = rest_I[k] + KYBER_N;
                e_u_x[k] = rest_I[k] + KYBER_N;
                s_u_y[k] = pi->u_s_2ddeg_shares[k][i][j];
                e_u_y[k] = pi->u_e_2ddeg_shares[k][i][j];
            }

            Vec<ZZ_p> ntl_s_u_x, ntl_s_u_y, ntl_e_u_x, ntl_e_u_y;
            ntl_s_u_x.SetLength(DEG_2D + 1);
            ntl_s_u_y.SetLength(DEG_2D + 1);
            ntl_e_u_x.SetLength(DEG_2D + 1);
            ntl_e_u_y.SetLength(DEG_2D + 1);
            for (size_t k = 0; k < DEG_2D + 1; ++k) {
                ntl_s_u_x[k] = s_u_x[k];
                ntl_e_u_x[k] = e_u_x[k];

                ntl_s_u_y[k] = s_u_y[k];
                ntl_e_u_y[k] = e_u_y[k];
            }
            ZZ_pX s_u_poly, e_u_poly;
            interpolate(s_u_poly, ntl_s_u_x, ntl_s_u_y);
            interpolate(e_u_poly, ntl_e_u_x, ntl_e_u_y);
            // eval the first KYBER_N points for eta shares for s / e:
            Vec<ZZ_p> s_u_y_ntl, e_u_y_ntl;
            e_u_y_ntl.SetLength(KYBER_N), s_u_y_ntl.SetLength(KYBER_N);
            for (size_t k = 0; k < KYBER_N; ++k) {
                s_u_y_ntl[k] = eval(s_u_poly, ZZ_p(k));
                e_u_y_ntl[k] = eval(e_u_poly, ZZ_p(k));
                // check u = 0 by interpolation points from N\I
                if (s_u_y_ntl[k] != gf3329_from_int(0)) {
                    printf("Check failed for s.u[%zu] at %zu (%d != 0).\n",
                           i, k, conv<int>(s_u_y_ntl[k]));
                    return false;
                }
                if (e_u_y_ntl[k] != gf3329_from_int(0)) {
                    printf("Check failed for e.u[%zu] at %zu (%d != 0).\n",
                           i, k, conv<int>(e_u_y_ntl[k]));
                    return false;
                }

            }

            // rebuild the shares of u_2d for s and e
            // if not in I, then set the value to as in pi.u_s_2ddeg_shares / u_e_2ddeg_shares
            for (size_t k = 0; k < MPCITH_N - MPCITH_T; ++k) {
                s_u_shares[i][j].share_x[rest_I[k]] = rest_I[k] + KYBER_N;
                s_u_shares[i][j].share_y[rest_I[k]] = pi->u_s_2ddeg_shares[k][i][j];
                e_u_shares[i][j].share_x[rest_I[k]] = rest_I[k] + KYBER_N;
                e_u_shares[i][j].share_y[rest_I[k]] = pi->u_e_2ddeg_shares[k][i][j];
            }
            // recon the 2d-degree shares by combining the points from N\I and re-computed u_2d points in I
            secret_vec recom_zero_s, recom_zero_e;
            recon_secrets_2ddeg(&recom_zero_s, &s_u_shares[i][j]);
            recon_secrets_2ddeg(&recom_zero_e, &e_u_shares[i][j]);
            // check the recom_zero_s, recom_zero_e are zero vectors
            for (size_t k = 0; k < KYBER_N; ++k) {
                if (recom_zero_s.secret[k] != 0) {
                    printf("Check failed for inconsistency of s.u2d[%zu] shares at %zu.\n",
                           i, k);
                    return false;
                }
                if (recom_zero_e.secret[k] != 0) {
                    printf("Check failed for inconsistency of e.u2d[%zu] shares at %zu.\n",
                           i, k);
                    return false;
                }
            }
        }
    }
    // set u_2d shares for parties in I
    uint16_t s_u2d_shares_view[MPCITH_T][KYBER_K][2 * KYBER_ETA1], e_u2d_shares_view[MPCITH_T][KYBER_K][2 * KYBER_ETA1];
    for (size_t i = 0; i < KYBER_K; ++i) {
        for (size_t j = 0; j < 2 * KYBER_ETA1; ++j) {
            for (size_t k = 0; k < MPCITH_T; ++k) {
                s_u2d_shares_view[k][i][j] = s_u_shares[i][j].share_y[pi->I[k]];
                e_u2d_shares_view[k][i][j] = e_u_shares[i][j].share_y[pi->I[k]];
            }
        }
    }

    // validate the commitments of comm.
    uint8_t ch_seeds[MPCITH_T][KYBER_SYMBYTES];
    for (size_t i = 0; i < MPCITH_T; ++i) {
        uint8_t view[KYBER_SYMBYTES + (((6 + 4 * KYBER_ETA1 * 2) * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1))) * sizeof(uint16_t)] = {0};
        // copy the views in set I:
        memcpy(view,
               Tcomm_rec[pi->I[i]],
               KYBER_SYMBYTES);
        memcpy(view + KYBER_SYMBYTES,
               pi->s_shares[i],
               KYBER_K * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + KYBER_K * sizeof(uint16_t),
               pi->e_shares[i],
               KYBER_K * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + (2 * KYBER_K) * sizeof(uint16_t),
               pi->f_shares[i],
               (MPCITH_K + MPCITH_V + 1) * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + (2 * KYBER_K + (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t),
               pi->NTT_f_shares[i],
               (MPCITH_K + MPCITH_V + 1) * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + (2 * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t),
               beta[pi->I[i]],
               KYBER_K * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + (3 * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t),
               gamma[pi->I[i]],
               KYBER_K * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + (4 * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t),
               sr_shares_view[i],
               KYBER_K * sizeof(uint16_t));
        memcpy(view + KYBER_SYMBYTES + (5 * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1)) * sizeof(uint16_t),
               er_shares_view[i],
               KYBER_K * sizeof(uint16_t));

        for (size_t j = 0; j < KYBER_K; ++j) {
            size_t curr_offset = KYBER_SYMBYTES + (6 * KYBER_K + 2 * (MPCITH_K + MPCITH_V + 1) + 4 * (j * 2 * KYBER_ETA1)) * sizeof(uint16_t);
            memcpy(view + curr_offset,
                   pi->z_s_ddeg_shares[i][j],
                   2 * KYBER_ETA1 * sizeof(uint16_t));
            memcpy(view + curr_offset + 1 * (2 * KYBER_ETA1) * sizeof(uint16_t),
                   pi->z_e_ddeg_shares[i][j],
                   2 * KYBER_ETA1 * sizeof(uint16_t));
            memcpy(view + curr_offset + 2 * (2 * KYBER_ETA1) * sizeof(uint16_t),
                   s_u2d_shares_view[i][j],
                   2 * KYBER_ETA1 * sizeof(uint16_t));
            memcpy(view + curr_offset + 3 * (2 * KYBER_ETA1) * sizeof(uint16_t),
                   e_u2d_shares_view[i][j],
                   2 * KYBER_ETA1 * sizeof(uint16_t));
        }
        sha3_256(ch_seeds[i], view, sizeof(view));
    }
    uint8_t merged_ch[MPCITH_N * KYBER_SYMBYTES], ch[KYBER_SYMBYTES];
    for (size_t i = 0, in_set = 0, out_set = 0; i < MPCITH_N; ++i) {
        if (in_I[i]) {
            size_t pos = 0;
            for (size_t j = 0; j < MPCITH_T; ++j) {
                if (pi->I[j] == i) {
                    pos = j;
                    break;
                }
            }
            memcpy(merged_ch + i * KYBER_SYMBYTES, ch_seeds[pos], KYBER_SYMBYTES);
            in_set += 1;
        } else {
            memcpy(merged_ch + i * KYBER_SYMBYTES, pi->comm[out_set], KYBER_SYMBYTES);
            out_set += 1;
        }
    }
    sha3_256(ch, merged_ch, MPCITH_N * KYBER_SYMBYTES);

    // recompute the index list of opened parties
    uint8_t I_[2 * MPCITH_T];
    uint16_t reom_I[MPCITH_T];
    kyber_shake256_prf(I_, 2 * MPCITH_T, ch, 1);
    // merge I_ into I and mod 1454
    for (size_t i = 0; i < MPCITH_T; ++i) {
        reom_I[i] = ((I_[2 * i] << 8) | I_[2 * i + 1]) % MPCITH_N;
    }
    // remove duplicate
    for (size_t i = 1; i < MPCITH_T; i++) {
        size_t j = 0;
        uint16_t inc = 0;
        bool is_dup;
        do {
            is_dup = false;
            for (j = 0; j < i; j++) {
                if ((reom_I[i] + inc) % MPCITH_N == reom_I[j]) {
                    is_dup = true;
                    inc = inc + 1;
                    break;
                }
            }
        } while (is_dup);
        reom_I[i] = (reom_I[i] + inc) % MPCITH_N;
    }
    // validate the reom_I is equal to pi->I
    for (size_t i = 0; i < MPCITH_T; ++i) {
        if (reom_I[i] != pi->I[i]) {
            printf("Check failed for reom_I[%zu]=%d (pi.I[%zu]=%d).\n", i, reom_I[i], i, pi->I[i]);
            return false;
        }
    }

    return true;
}