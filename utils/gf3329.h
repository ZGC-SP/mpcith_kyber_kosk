
#ifndef MPCITH_KYBER_GF3329_H
#define MPCITH_KYBER_GF3329_H

#define MODULUS 3329

#include <stdint.h>

#include "../kyber/reduce.h"

uint16_t gf3329_add(uint16_t a, uint16_t b); // return a+b
uint16_t gf3329_sub(uint16_t a, uint16_t b); // return a-b
uint16_t gf3329_mul(uint16_t a, uint16_t b); // return a*b
uint16_t gf3329_pow(uint16_t a, uint16_t b); // return a^b
uint16_t gf3329_inv(uint16_t a); // return a^(-1)

// encode from int16_t to uint16_t (gf3329)
uint16_t encode_to_gf3329(int16_t a);

// decode from uint16_t (gf3329) to int16_t
int16_t decode_from_gf3329(uint16_t a);

// encode from int to uint16_t (gf3329)
uint16_t encode_to_gf3329_2(int a);

// decode from uint16_t (gf3329) to int
int decode_from_gf3329_2(uint16_t a);

int16_t gf3329_int_to_int16(int a);

uint16_t gf3329_from_int(int a);


#endif //MPCITH_KYBER_GF3329_H
