//
// Created by 高谦文 on 2023/4/21.
//

#ifndef PSS_FRODO_PROVER_H
#define PSS_FRODO_PROVER_H

#include <NTL/ZZ_p.h>
#include <NTL/matrix.h>
#include <NTL/pair.h>

#include "packed_secret_sharing.h"
#include "plain_LWE_prover.h"
#include "utils.h"
#include <thread>

class Frodo_Prover {
    /**
     * Prove the relation AS + E = B in Frodo KEM.
     */

private:
    long n{}; // the size of A (n * n)
    long m{}; // the size of S & E (n * m)

    long k{}; // the security parameter (for the soundness error)
    long v{}; // the batch size 8

    Mat<ZZ_p> A; // the matrix A

    Vec<ZZ_p> R;
    Vec<ZZ_p> AR;

    Vec<Pair<ZZ_p, ZZ_p>> R_shares; // the random vector r for pre-processing
    Vec<Pair<ZZ_p, ZZ_p>> AR_shares; // the A * r for pre-processing


    Vec<Vec<Pair<ZZ_p, ZZ_p>>> const_share_s, const_share_e;
    PSS::PackedSecretSharing pss; // the packed secret-sharing parameters

    long default_thread_num;



public:

    Frodo_Prover(long n, long m,
                 long k,
                 long t, long N,
                 const Mat<ZZ_p> &mat_A,
                 long default_thread_num = 1);

    ~Frodo_Prover() = default;

    void InitializeConstantShares();
    // unit prove: prove the relation: A * s + e = b; s, e, b are one row/col in S, E, B.
    static void UnitProve(const Frodo_Prover& frodo, long index, const Vec<ZZ_p> &s, const Vec<ZZ_p> &e, const Vec<ZZ_p> &b);
    static void Prove(const Frodo_Prover& frodo, const Mat<ZZ_p> &mat_S, const Mat<ZZ_p> &mat_E, const Mat<ZZ_p> &mat_B);



};


#endif //PSS_FRODO_PROVER_H
