
#ifndef MPCITH_KYBER_KOSK_MLWE_VERIFIER_HPP
#define MPCITH_KYBER_KOSK_MLWE_VERIFIER_HPP

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/BasicThreadPool.h>
using namespace NTL;

#include "params.hpp"
#include "mlwe_prover.hpp"

bool verify(const mpcith_proof *pi,
            const mlwe_inst *mlwe);

#endif //MPCITH_KYBER_KOSK_MLWE_VERIFIER_HPP
