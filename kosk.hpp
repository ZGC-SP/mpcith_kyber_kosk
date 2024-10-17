#ifndef MPCITH_KYBER_KOSK_KOSK_HPP
#define MPCITH_KYBER_KOSK_KOSK_HPP

extern "C" {
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
}

#include "mlwe_prover.hpp"
#include "mlwe_verifier.hpp"

typedef struct {
    uint8_t pk[KYBER_PUBLICKEYBYTES];
    uint8_t sk[KYBER_SECRETKEYBYTES];
} kyber_keypair;

void kyber_keygen(kyber_keypair *keypair,
                  mlwe_inst *raw_key);
void kyber_verifiable_keygen(kyber_keypair *keypair,
                             uint8_t *pi);

bool kyber_kosk_verify(const uint8_t *pi,
                       const uint8_t *pk);



#endif //MPCITH_KYBER_KOSK_KOSK_HPP
