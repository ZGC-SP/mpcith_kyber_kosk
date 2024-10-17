#ifndef MPCITH_KYBER_KOSK_SS_HPP
#define MPCITH_KYBER_KOSK_SS_HPP

extern "C" {
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include "kyber/randombytes.h"
#include "utils/gf3329.h"
#include "utils/precomputed_kyber.h"
}

#include "params.hpp"

/**
 * the share polynomial format:
 * coefficients:
 * -----------------------------------------------------------------------------
 * |f[0] f[1] f[2] ... f[kyber_N-1]| f[kyber_N] ... f[kyber_N+t] f[kyber_N+t+1]| // degree = (256 + t)
 * |---------------------------------------------------------------------------|
 * | embed kyber-kem secret (256)  |               random (t+1)                |
 * -----------------------------------------------------------------------------
 *
 * shares:
 * --------------------------------------------------------------------------------------------------
 * |F(0) F(1) F(2) ... F(kyber_N-1)| F(kyber_N) ... F(kyber_N+t) F(kyber_N+t+1) ... F(mpc_N+kyber_N)|
 * |------------------------------------------------------------------------------------------------|
 * | embed kyber-kem secret (256)  |                 shares for MPCITH_N parties                    |
 * |     not share this offset     |                the x-coord start from kyber_N                  |
 * --------------------------------------------------------------------------------------------------
 */

typedef struct {
    size_t len;
    uint16_t share_x[MPCITH_N];
    uint16_t share_y[MPCITH_N];
} share_vec;

typedef struct {
    size_t len;
    uint16_t secret[KYBER_N];
} secret_vec;

int share_secrets_ddeg(share_vec *shares, const secret_vec *secrets);

int recon_secrets_ddeg(secret_vec *secrets, const share_vec *shares);
int recon_secrets_2ddeg(secret_vec *secrets, const share_vec *shares);

int recompute_share_secrets_ddeg(share_vec *shares, const uint16_t *y);
//int recompute_share_secrets_2ddeg(share_vec *shares, const uint16_t *y);

int shares_add(share_vec *res, const share_vec *a, const share_vec *b);
int shares_sub(share_vec *res, const share_vec *a, const share_vec *b);
int shares_mul(share_vec *res_2ddeg, const share_vec *a, const share_vec *b);

#define DEG_D (KYBER_N + MPCITH_T)
#define DEG_2D ((KYBER_N + MPCITH_T) * 2)
#define EVAL_POINTS KYBER_N

void gen_pow_table(uint16_t *tab, uint16_t x, uint16_t deg);
void interpolate(uint16_t* coeffs, const uint16_t* x, const uint16_t* y, size_t n);
uint16_t poly_eval(const uint16_t *coeff, const uint16_t *x_pow_table, uint16_t deg);
void poly_eval_from_points_ddeg(uint16_t *res, const uint16_t *x, const uint16_t *y);
void poly_eval_from_points_2ddeg(uint16_t *res, const uint16_t *x, const uint16_t *y);
void poly_eval_from_points(uint16_t *res,
                           const uint16_t *x, const uint16_t *y, uint16_t deg,
                           const uint16_t *eval_x, uint16_t eval_n);

#endif //MPCITH_KYBER_KOSK_SS_HPP
