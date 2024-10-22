//
// Created by 高谦文 on 2023/4/27.
//

#ifndef PSS_KYBER_PROVER_H
#define PSS_KYBER_PROVER_H

#include <NTL/ZZ_p.h>
#include <NTL/pair.h>
#include <NTL/matrix.h>

#include "packed_secret_sharing.h"
#include "plain_LWE_prover.h"
#include "utils.h"

class KyberProver {

private:
    long n{};
    long k{};
    long v{}; // v = 4

    Mat<ZZ_p> T; // the NTT matrix

    Vec<ZZ_p> R;
    Vec<ZZ_p> TR;

    Vec<Pair<ZZ_p, ZZ_p>> R_shares;
    Vec<Pair<ZZ_p, ZZ_p>> TR_shares;

    Vec<Vec<Pair<ZZ_p, ZZ_p>>> const_share_s;
    Vec<Vec<Pair<ZZ_p, ZZ_p>>> const_share_e;

    long default_thread_num{};

    PSS::PackedSecretSharing pss;
public:
    KyberProver(long n, long N, long t, long k, const Mat<ZZ_p> &mat_T, long default_thread_num = 1);

    ~KyberProver() = default;

    void InitializeConstantShares();

    static void Prove(const KyberProver &kyber,
                      const Mat<ZZ_pX> &A, const Mat<ZZ_p> &S, const Mat<ZZ_p> &E, const Mat<ZZ_p> &B);

    static Vec<Vec<Pair<ZZ_p, ZZ_p>>> ShortVectorProve(const KyberProver &kyber,
                                 const Vec<Pair<ZZ_p, ZZ_p>> &share_s, const Vec<Pair<ZZ_p, ZZ_p>> &share_e);
};


#endif //PSS_KYBER_PROVER_H
