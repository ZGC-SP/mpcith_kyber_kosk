#include "ss.hpp"

int share_secrets_ddeg(share_vec *shares, const secret_vec *secrets) {
    uint8_t random_bytes[(MPCITH_T + 1) * 2];
    randombytes(random_bytes, (MPCITH_T + 1) * 2);
    uint16_t random_uint16[MPCITH_T + 1];
    for (size_t i = 0; i < MPCITH_T + 1; ++i) {
        random_uint16[i] = ((random_bytes[i * 2] << 8) | random_bytes[i * 2 + 1]) % KYBER_Q;
        shares->share_x[i] = i + KYBER_N;
        shares->share_y[i] = random_uint16[i];
    }

    uint16_t sec[KYBER_N + MPCITH_T + 1] = {0};
    for (size_t i = 0; i < KYBER_N + MPCITH_T + 1; ++i) {
        if (i < KYBER_N) {
            sec[i] = secrets->secret[i];
        } else {
            sec[i] = random_uint16[i - KYBER_N];
        }
    }

    // evaluate the rest shares
    for (size_t i = MPCITH_T + KYBER_N + 1; i < MPCITH_N + KYBER_N; ++i) {
        uint16_t share_y = 0;
        for (size_t j = 0; j < MPCITH_T + KYBER_N + 1; ++j) {
            uint16_t prod = gf3329_mul(sec[j],
                                       get_precomputed_share_coeff_ddeg(i - (MPCITH_T + KYBER_N + 1), j));
            share_y = gf3329_add(share_y, prod);
        }
        shares->share_x[i - KYBER_N] = i;
        shares->share_y[i - KYBER_N] = share_y;
    }
    return 0;
}


int recon_secrets_ddeg(secret_vec *secrets, const share_vec *shares) {
    int share_offset = KYBER_N;
    int degree = KYBER_N + MPCITH_T;
    if (shares->share_x[0] != share_offset) {
        // error
        return -1;
    }
    for (size_t i = 0; i < KYBER_N; ++i) {
        uint16_t rec_s = 0;
        for (size_t j = 0; j < degree + 1; ++j) {
            uint16_t prod = gf3329_mul(shares->share_y[j], get_precomputed_recon_coeff_ddeg(i, j));
            rec_s = gf3329_add(rec_s, prod);
        }
        secrets->secret[i] = rec_s;
    }

    return 0;
}

int recon_secrets_2ddeg(secret_vec *secrets, const share_vec *shares) {
    int share_offset = KYBER_N;
    int degree = 2 * (KYBER_N + MPCITH_T);
    if (shares->share_x[0] != share_offset) {
        // error
        return -1;
    }
    for (int i = 0; i < KYBER_N; ++i) {
        uint16_t rec_s = 0;
        for (int j = 0; j < degree + 1; ++j) {
            uint16_t prod = gf3329_mul(shares->share_y[j], get_precomputed_recon_coeff_2ddeg(i, j));
            rec_s = gf3329_add(rec_s, prod);
        }
        secrets->secret[i] = rec_s;
    }

    return 0;
}

// here the length of y vec is KYBER_N + MPCITH_T + 1 (i.e. degree of the polynomial)
int recompute_share_secrets_ddeg(share_vec *shares, const uint16_t *y) {
    for (size_t i = 0; i < MPCITH_T + 1; ++i) {
        shares->share_x[i] = i + KYBER_N;
        shares->share_y[i] = y[i + KYBER_N];
    }

    uint16_t sec[KYBER_N + MPCITH_T + 1] = {0};
    for (size_t i = 0; i < KYBER_N + MPCITH_T + 1; ++i) {
        sec[i] = y[i];
    }

    // evaluate the rest shares
    for (size_t i = MPCITH_T + KYBER_N + 1; i < MPCITH_N + KYBER_N; ++i) {
        uint16_t share_y = 0;
        for (size_t j = 0; j < MPCITH_T + KYBER_N + 1; ++j) {
            uint16_t prod = gf3329_mul(sec[j],
                                       get_precomputed_share_coeff_ddeg(i - (MPCITH_T + KYBER_N + 1), j));
            share_y = gf3329_add(share_y, prod);
        }
        shares->share_x[i - KYBER_N] = i;
        shares->share_y[i - KYBER_N] = share_y;
    }
    return 0;
}

int shares_add(share_vec *res, const share_vec *a, const share_vec *b) {
    for (size_t i = 0; i < MPCITH_N; ++i) {
        if (a->share_x[i] != b->share_x[i]) {
            // error
            return -1;
        }
        res->share_x[i] = a->share_x[i];
        res->share_y[i] = gf3329_add(a->share_y[i], b->share_y[i]);
    }
    return 0;
}


int shares_sub(share_vec *res, const share_vec *a, const share_vec *b) {
    for (size_t i = 0; i < MPCITH_N; ++i) {
        if (a->share_x[i] != b->share_x[i]) {
            // error
            return -1;
        }
        res->share_x[i] = a->share_x[i];
        res->share_y[i] = gf3329_sub(a->share_y[i], b->share_y[i]);
    }
    return 0;
}

int shares_mul(share_vec *res_2ddeg, const share_vec *a, const share_vec *b) {
    for (size_t i = 0; i < MPCITH_N; ++i) {
        if (a->share_x[i] != b->share_x[i]) {
            // error
            return -1;
        }
        res_2ddeg->share_x[i] = a->share_x[i];
        res_2ddeg->share_y[i] = gf3329_mul(a->share_y[i], b->share_y[i]);
    }
    return 0;
}

void interpolate(uint16_t* coeffs, const uint16_t* x, const uint16_t* y, size_t n) {
    for (int i = 0; i < n; i++) {
        coeffs[i] = 0; // Initialize all coefficients to 0
    }

    // Iterate through each point (xi, yi)
    for (int i = 0; i < n; i++) {
        uint16_t numerator[n], denominator = 1;

        // Initialize numerator to the constant term 1
        numerator[0] = 1;
        for (int j = 1; j < n; j++) {
            numerator[j] = 0;
        }

        // Compute the Lagrange basis polynomial L_i(x)
        for (int j = 0; j < n; j++) {
            if (i != j) {
                denominator = gf3329_mul(denominator,
                                         gf3329_add(x[i],
                                                    gf3329_sub(MODULUS, x[j])));

                // Multiply the numerator by (x - xj)
                for (int k = n - 1; k >= 0; k--) {
                    numerator[k] = gf3329_mul(numerator[k], gf3329_sub(MODULUS, x[j]));
                    if (k > 0) {
                        numerator[k] = gf3329_add(numerator[k], numerator[k - 1]);
                    }
                }
            }
        }

        // Invert the denominator
        denominator = gf3329_inv(denominator);

        // Multiply the polynomial by yi / denominator
        for (int k = 0; k < n; k++) {
            coeffs[k] = gf3329_add(coeffs[k], gf3329_mul(numerator[k], gf3329_mul(y[i], denominator)));
        }
    }
}

void gen_pow_table(uint16_t *tab, uint16_t x, uint16_t deg) {
    tab[0] = 1;
    tab[1] = x;
    for (size_t i = 2; i < deg + 1; ++i ) {
        tab[i] = gf3329_mul(tab[i - 1], x);
    }
}

uint16_t poly_eval(const uint16_t *coeff, const uint16_t *x_pow_table, uint16_t deg) {
    uint16_t r = 0;
    for (size_t i = 0; i < deg + 1; ++i ) {
        r = gf3329_add(r, gf3329_mul(coeff[i], x_pow_table[i]));
    }
    return r;
}

void poly_eval_from_points_ddeg(uint16_t *res, const uint16_t *x, const uint16_t *y) {
    bool already_in[DEG_D + 1] = {false};
    for (int i = 0; i < DEG_D + 1; i++) {
        if (x[i] < DEG_D + 1) {
            already_in[x[i]] = true;
        }
    }

    for (int i = 0; i < DEG_D + 1; i++) {
        if (already_in[i]) {
            res[i] = y[i];
        } else {
            uint16_t result = 0;
            for (int j = 0; j < DEG_D + 1; j++) {
                uint16_t numerator = 1;
                uint16_t denominator = 1;
                for (int k = 0; k < DEG_D + 1; k++) {
                    if (k != j) {
                        numerator = gf3329_mul(numerator, gf3329_sub(i, x[k]));
                        denominator = gf3329_mul(denominator, gf3329_sub(x[j], x[k]));
                    }
                }
                uint16_t term = gf3329_mul(y[j], gf3329_mul(numerator, gf3329_inv(denominator)));
                result = gf3329_add(result, term);
            }
            res[i] = result;
        }
    }
}

void poly_eval_from_points_2ddeg(uint16_t *res, const uint16_t *x, const uint16_t *y) {
    for (int i = 0; i < DEG_2D + 1; i++) {
        uint16_t result = 0;
        for (int j = 0; j < DEG_2D + 1; j++) {
            uint16_t numerator = 1;
            uint16_t denominator = 1;
            for (int k = 0; k < DEG_2D + 1; k++) {
                if (k != j) {
                    numerator = gf3329_mul(numerator, gf3329_sub(i, x[k]));
                    denominator = gf3329_mul(denominator, gf3329_sub(x[j], x[k]));
                }
            }
            uint16_t term = gf3329_mul(y[j], gf3329_mul(numerator, gf3329_inv(denominator)));
            result = gf3329_add(result, term);
        }
        res[i] = result;
    }
}

void poly_eval_from_points(uint16_t *res,
                           const uint16_t *x, const uint16_t *y, uint16_t deg,
                           const uint16_t *eval_x, uint16_t eval_n) {
    for (int i = 0; i < eval_n; i++) {
        uint16_t result = 0;
        for (int j = 0; j < deg + 1; j++) {
            uint16_t numerator = 1;
            uint16_t denominator = 1;
            for (int k = 0; k < deg + 1; k++) {
                if (k != j) {
                    numerator = gf3329_mul(numerator, gf3329_sub(eval_x[i], x[k]));
                    denominator = gf3329_mul(denominator, gf3329_sub(x[j], x[k]));
                }
            }
            uint16_t term = gf3329_mul(y[j], gf3329_mul(numerator, gf3329_inv(denominator)));
            result = gf3329_add(result, term);
        }
        res[i] = result;
    }
}