//
// Created by 高谦文 on 2024/7/8.
//

#ifndef MPCITH_KYBER_PRECOMPUTED_KYBER512_H
#define MPCITH_KYBER_PRECOMPUTED_KYBER512_H

#include <stdint.h>

uint16_t get_precomputed_share_coeff_ddeg(int x, int i);
uint16_t get_precomputed_share_coeff_2ddeg(int x, int i);
uint16_t get_precomputed_recon_coeff_ddeg(int x, int i);
uint16_t get_precomputed_recon_coeff_2ddeg(int x, int i);

#endif //MPCITH_KYBER_PRECOMPUTED_KYBER512_H
