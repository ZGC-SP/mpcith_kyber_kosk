//
// Created by 高谦文 on 2023/4/14.
//

#ifndef PSS_PLAIN_LWE_PROVER_H
#define PSS_PLAIN_LWE_PROVER_H

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/pair.h>
#include <vector>
#include <NTL/vector.h>
#include "packed_secret_sharing.h"
#include "utils.h"

using namespace NTL;

struct PlainLWE {
    Mat<ZZ_p> A;
    std::vector<Vec<ZZ_p>> s, e, b;
};

PlainLWE GeneratePlainLweInstance(long m, long n, long size);

class PlainLWE_Prover {
    /**
     * Prove the relation: As + e = b, where:
     *  A is a m*n matrix, \in Z_p
     *  s is a n*1 vector, \in Z_2,
     *  b is a m*1 vector, \in Z_p,
     *  e is a m*1 vector(noise), \in Z_2
     */

private:
    // internal variables

    long m{}; // the row number of A
    long n{}; // the column number of A (also the secret number of PSS)
    long k{}; // the security parameter (for the soundness error)
    long v{}; // the batch size

    Mat<ZZ_p> A; // the matrix A
    Vec<ZZ_p> r; // the random vector r for pre-processing
    Vec<ZZ_p> Ar; // the A * r for pre-processing

    Vec<Vec<Pair<ZZ_p, ZZ_p>>> vec_f_shares, vec_Af_shares; // the shares of each f and Af (k+v items)

    // the random vector r is evaluate by: F(a) = a^0 * vec_f[0] + a^1 * vec_f[1] + ... + a^(k+v) * vec_f[k+v]
    Vec<Pair<ZZ_p, ZZ_p>> r_shares, Ar_shares; // the shares of r and Ar

    // the Ar is evaluate by: F(a) = a^0 * vec_Af[0] + a^1 * vec_Af[1] + ... + a^(k+v) * vec_Af[k+v]
    // the shares of all-1 vector [1]
    Vec<Pair<ZZ_p, ZZ_p>> const_share_s, const_share_e;

    // The commitment of the shares of f and A * f
    Mat<std::string> commitments;

    PSS::PackedSecretSharing pss; // the packed secret-sharing parameters

    long default_thread_num;

public:
    // constructor & destructor
    /**
     * Initialize the parameters of the the plain LWE prover (A*s + e = b)
     * @param k the security parameter (for the soundness error)
     * @param v the batch size
     * @param t the threshold of adversary (for MPC-it-Head)
     * @param N the number of parties for packed-secret sharing
     * @param mat_A the Matrix of LWE.
     */
    PlainLWE_Prover(long k, long v, long t, long N, const Mat<ZZ_p> &mat_A, long default_thread_num = 1);

    ~PlainLWE_Prover();

    // getter
    long GetM() const { return this->m; }
    long GetN() const { return this->n; }
    long GetK() const { return this->k; }
    long GetV() const { return this->v; }
    Vec<ZZ_p> GetR() const { return this->r; }
    Vec<ZZ_p> GetAr() const { return this->Ar; }

    Mat<ZZ_p> GetA() const { return this->A; }
    Vec<Pair<ZZ_p, ZZ_p>> GetRShares() const { return this->r_shares; }
    Vec<Pair<ZZ_p, ZZ_p>> GetArShares() const { return this->Ar_shares; }

    // generate the constant shares to prove s and e are short
    // generate two shares of vector [1]
    void InitializeConstantShares();

    void Prove(Vec<ZZ_p> s, Vec<ZZ_p> e, Vec<ZZ_p> b);

};


#endif //PSS_PLAIN_LWE_PROVER_H
