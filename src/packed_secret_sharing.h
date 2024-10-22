//
// Created by 高谦文 on 2023/4/6.
//

#ifndef PSS_PACKED_SECRET_SHARING_H
#define PSS_PACKED_SECRET_SHARING_H

#include <vector>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>

#include "utils.h"

namespace PSS {

    // FIXME: use template to support different domain of secret (over Z_q and GF(2^m))
    class PackedSecretSharing {
    private:
        // internal variables
        ZZ p; // prime number of field Z_p
        long N{}; // number of parties
        long t{}; // adversary threshold
        long l{}; // number of shares
        long d{}; // degree of the secret polynomial

        // interpolation optimized table
        Vec<ZZ_p> factorial_table;
        Vec<ZZ_pX> lagrange_table;

    public:
        // PSS constructor & destructor
        PackedSecretSharing();

        PackedSecretSharing(long N, long t, long l);

        ~PackedSecretSharing();

        // getters & setters
        long GetN() const { return this->N; }
        long GetThreshold() const { return this->t; }
        long GetSecretNumber() const { return this->l; }
        long GetDegree() const { return this->d; }

        /// Core implementations
        // PSS methods
        Vec<Pair<ZZ_p, ZZ_p>> Share(
                const Vec<ZZ_p> &secrets,
                long thread_num = 1,
                bool isContinuous = false) const;

        Vec<ZZ_p> Reconstruct(
                const Vec<Pair<ZZ_p, ZZ_p>> &shares,
                long thread_num = 1,
                bool isContinuous = false) const;

        // FIXME: fast share & reconstruct, confirmation needed
        Vec<Pair<ZZ_p, ZZ_p>> FastShare(
                const Vec<ZZ_p> &secrets,
                long thread_num = 1) const;
        Vec<ZZ_p> FastReconstruct(
                const Vec<Pair<ZZ_p, ZZ_p>> &shares,
                long thread_num = 1) const;

        // Evaluations
        // z = x + y
        Vec<Pair<ZZ_p, ZZ_p>> Addition(
                const Vec<Pair<ZZ_p, ZZ_p>> &x,
                const Vec<Pair<ZZ_p, ZZ_p>> &y,
                long thread_num = 1) const;

        // z = x - y
        Vec<Pair<ZZ_p, ZZ_p>> Subtraction(
                const Vec<Pair<ZZ_p, ZZ_p>> &x,
                const Vec<Pair<ZZ_p, ZZ_p>> &y,
                long thread_num = 1) const;

        // z = x * y
        Vec<Pair<ZZ_p, ZZ_p>> Multiplication(
                const Vec<Pair<ZZ_p, ZZ_p>> &x,
                const Vec<Pair<ZZ_p, ZZ_p>> &y,
                long thread_num = 1) const;

        Vec<Pair<ZZ_p, ZZ_p>> MultiplicationReShare(
                const Vec<Pair<ZZ_p, ZZ_p>> &mul_shares,
                long thread_num = 1,
                bool isContinuous = false) const;
    };

} // PSS

#endif //PSS_PACKED_SECRET_SHARING_H
