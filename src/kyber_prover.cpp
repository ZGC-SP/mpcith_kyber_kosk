//
// Created by 高谦文 on 2023/4/27.
//

#include "kyber_prover.h"

KyberProver::KyberProver(long n, long N, long t, long k, const Mat<ZZ_p> &T, long default_thread_num) {
    this->T = T;
    this->n = n;
    this->k = k;
    this->v = 4;

    this->pss = PSS::PackedSecretSharing(N, t, this->n);

    this->default_thread_num = default_thread_num;

    PlainLWE_Prover pLWE = PlainLWE_Prover(this->k, this->v, t, N, this->T, default_thread_num);

    this->R = pLWE.GetR();
    this->TR = pLWE.GetAr();

    this->R_shares = pLWE.GetRShares();
    this->TR_shares = pLWE.GetArShares();

}

void KyberProver::InitializeConstantShares() {
    // the element value of s and e is ranged from [0, +/-12]
    long range = 2;
    this->const_share_s.SetLength(2 * range + 1);
    this->const_share_e.SetLength(2 * range + 1);
    Vec<Vec<ZZ_p>> constant_value;
    constant_value.SetLength(2 * range + 1);
    for (int i = 0; i < 2 * range + 1; ++i) {
        constant_value[i].SetLength(this->n, ZZ_p(i));
        Vec<ZZ_p> range_vec;
        range_vec.SetLength(this->n, ZZ_p(range));
        Vec<ZZ_p> vec_to_share = constant_value[i] - range_vec;

        this->const_share_s[i] = this->pss.Share(vec_to_share);
        this->const_share_e[i] = this->pss.Share(vec_to_share);
    }
}

void KyberProver::Prove(const KyberProver &kyber,
                        const Mat<ZZ_pX> &A, const Mat<ZZ_p> &S, const Mat<ZZ_p> &E, const Mat<ZZ_p> &B) {
    /**
     * A = [a00 a01]  , S = [s0] , E = [e0]
     *     [a10 a11]        [s1]       [e1]
     */
    long range = 2;

    auto start_linear_proof = std::chrono::high_resolution_clock::now();
    // share s0, s1, e0, e1:
    Vec<Pair<ZZ_p, ZZ_p>> shares_s0 = kyber.pss.Share(S[0], kyber.default_thread_num, true);
    Vec<Pair<ZZ_p, ZZ_p>> shares_s1 = kyber.pss.Share(S[1], kyber.default_thread_num, true);

    Vec<Pair<ZZ_p, ZZ_p>> shares_e0 = kyber.pss.Share(E[0], kyber.default_thread_num, true);
    Vec<Pair<ZZ_p, ZZ_p>> shares_e1 = kyber.pss.Share(E[1], kyber.default_thread_num, true);

    // Phase 1: Prove the linear relation
    {
        // NTT(a[i,j])
        Vec<ZZ_p> NTT_a00 = kyber.T * PolynomialToVector(A.get(0, 0), kyber.n);
        Vec<ZZ_p> NTT_a01 = kyber.T * PolynomialToVector(A.get(0, 1), kyber.n);
        Vec<ZZ_p> NTT_a10 = kyber.T * PolynomialToVector(A.get(1, 0), kyber.n);
        Vec<ZZ_p> NTT_a11 = kyber.T * PolynomialToVector(A.get(1, 1), kyber.n);

        Vec<Pair<ZZ_p, ZZ_p>> share_a00 = kyber.pss.Share(NTT_a00, kyber.default_thread_num, true);
        Vec<Pair<ZZ_p, ZZ_p>> share_a01 = kyber.pss.Share(NTT_a01, kyber.default_thread_num, true);
        Vec<Pair<ZZ_p, ZZ_p>> share_a10 = kyber.pss.Share(NTT_a10, kyber.default_thread_num, true);
        Vec<Pair<ZZ_p, ZZ_p>> share_a11 = kyber.pss.Share(NTT_a11, kyber.default_thread_num, true);

        // Handle:
        // 1. NTT(s0) * NTT(a00) + NTT(s1) * NTT(a01) + NTT(e0) = NTT(b0)
        // 2. NTT(s0) * NTT(a10) + NTT(s1) * NTT(a11) + NTT(e1) = NTT(b1)

        // handle s0, s1

        Vec<Pair<ZZ_p, ZZ_p>> share_s0_add_R = kyber.pss.Addition(shares_s0, kyber.R_shares, kyber.default_thread_num);
        Vec<Pair<ZZ_p, ZZ_p>> share_s1_add_R = kyber.pss.Addition(shares_s1, kyber.R_shares, kyber.default_thread_num);

        Vec<ZZ_p> vec_s0_add_R = kyber.pss.Reconstruct(share_s0_add_R, kyber.default_thread_num, true);
        Vec<ZZ_p> vec_s1_add_R = kyber.pss.Reconstruct(share_s1_add_R, kyber.default_thread_num, true);

        // NTT(s0 + R) - [R] , NTT(s1 + R) - [R]
        Vec<ZZ_p> NTT_s0_add_R = kyber.T * vec_s0_add_R;
        Vec<Pair<ZZ_p, ZZ_p>> share_NTT_s0_add_R = kyber.pss.Share(NTT_s0_add_R, kyber.default_thread_num, true);
        Vec<Pair<ZZ_p, ZZ_p>> share_NTT_s0 = kyber.pss.Subtraction(share_NTT_s0_add_R, kyber.TR_shares,
                                                                   kyber.default_thread_num);

        Vec<ZZ_p> NTT_s1_add_R = kyber.T * vec_s1_add_R;
        Vec<Pair<ZZ_p, ZZ_p>> share_NTT_s1_add_R = kyber.pss.Share(NTT_s1_add_R, kyber.default_thread_num, true);
        Vec<Pair<ZZ_p, ZZ_p>> share_NTT_s1 = kyber.pss.Subtraction(share_NTT_s1_add_R, kyber.TR_shares,
                                                                   kyber.default_thread_num);

        // handle e0, e1
        Vec<Pair<ZZ_p, ZZ_p>> share_e0_add_R = kyber.pss.Addition(shares_e0, kyber.R_shares, kyber.default_thread_num);
        Vec<Pair<ZZ_p, ZZ_p>> share_e1_add_R = kyber.pss.Addition(shares_e1, kyber.R_shares, kyber.default_thread_num);

        Vec<ZZ_p> vec_e0_add_R = kyber.pss.Reconstruct(share_e0_add_R, kyber.default_thread_num, true);
        Vec<ZZ_p> vec_e1_add_R = kyber.pss.Reconstruct(share_e1_add_R, kyber.default_thread_num, true);

        // NTT(s0 + R) - [R] , NTT(s1 + R) - [R]
        Vec<ZZ_p> NTT_e0_add_R = kyber.T * vec_e0_add_R;
        Vec<Pair<ZZ_p, ZZ_p>> share_NTT_e0_add_R = kyber.pss.Share(NTT_e0_add_R, kyber.default_thread_num, true);
        Vec<Pair<ZZ_p, ZZ_p>> share_NTT_e0 = kyber.pss.Subtraction(share_NTT_e0_add_R, kyber.TR_shares,
                                                                   kyber.default_thread_num);

        Vec<ZZ_p> NTT_e1_add_R = kyber.T * vec_e1_add_R;
        Vec<Pair<ZZ_p, ZZ_p>> share_NTT_e1_add_R = kyber.pss.Share(NTT_e1_add_R, kyber.default_thread_num, true);
        Vec<Pair<ZZ_p, ZZ_p>> share_NTT_e1 = kyber.pss.Subtraction(share_NTT_e1_add_R, kyber.TR_shares,
                                                                   kyber.default_thread_num);

        // build linear relation:
        Vec<Pair<ZZ_p, ZZ_p>> a00_mul_s0 = kyber.pss.Multiplication(share_a00, share_NTT_s0, kyber.default_thread_num);
        Vec<Pair<ZZ_p, ZZ_p>> a01_mul_s1 = kyber.pss.Multiplication(share_a01, share_NTT_s1, kyber.default_thread_num);

        Vec<Pair<ZZ_p, ZZ_p>> relation_0 = kyber.pss.Addition(a00_mul_s0, a01_mul_s1, kyber.default_thread_num);
        relation_0 = kyber.pss.Addition(relation_0, share_NTT_e0, kyber.default_thread_num);

        Vec<Pair<ZZ_p, ZZ_p>> a10_mul_s0 = kyber.pss.Multiplication(share_a10, share_NTT_s0, kyber.default_thread_num);
        Vec<Pair<ZZ_p, ZZ_p>> a11_mul_s1 = kyber.pss.Multiplication(share_a11, share_NTT_s1, kyber.default_thread_num);

        Vec<Pair<ZZ_p, ZZ_p>> relation_1 = kyber.pss.Addition(a10_mul_s0, a11_mul_s1, kyber.default_thread_num);

        relation_1 = kyber.pss.Addition(relation_1, share_NTT_e1, kyber.default_thread_num);
    }
    auto end_linear_proof = std::chrono::high_resolution_clock::now();

    // proof time:
    std::cout << "Short proof of linear relation time cost: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_linear_proof - start_linear_proof).count()
              << " ms" << std::endl;

    // Phase 2: prove the s and e are short:

    Vec<Vec<Pair<ZZ_p, ZZ_p>>> short_proof0;
    Vec<Vec<Pair<ZZ_p, ZZ_p>>> short_proof1;

    auto start_short_proof = std::chrono::high_resolution_clock::now();
    KyberProver::ShortVectorProve(kyber, shares_s0, shares_e0);
    KyberProver::ShortVectorProve(kyber, shares_s1, shares_e1);
    auto end_short_proof = std::chrono::high_resolution_clock::now();

    // proof time:
    std::cout << "Short proof of s & e time cost: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_short_proof - start_short_proof).count()
              << " ms" << std::endl;

    /**
     * return to: relation proof and short proof.
     */
}

Vec<Vec<Pair<ZZ_p, ZZ_p>>> KyberProver::ShortVectorProve(const KyberProver &kyber,
                                                         const Vec<Pair<ZZ_p, ZZ_p>> &share_s,
                                                         const Vec<Pair<ZZ_p, ZZ_p>> &share_e) {

    Vec<Pair<ZZ_p, ZZ_p>> share_U_s = kyber.pss.Subtraction(share_s, kyber.const_share_s[0], kyber.default_thread_num);
    Vec<Pair<ZZ_p, ZZ_p>> share_U_e = kyber.pss.Subtraction(share_e, kyber.const_share_e[0], kyber.default_thread_num);

    for (int i = 1; i < kyber.const_share_s.length(); ++i) {
        share_U_s = kyber.pss.Multiplication(share_U_s,
                                             kyber.pss.Subtraction(share_s, kyber.const_share_s[i],
                                                                   kyber.default_thread_num),
                                             kyber.default_thread_num);
        share_U_s = kyber.pss.MultiplicationReShare(share_U_s, kyber.default_thread_num);
    }
    for (int i = 1; i < kyber.const_share_e.length(); ++i) {
        share_U_e = kyber.pss.Multiplication(share_U_e,
                                             kyber.pss.Subtraction(share_e, kyber.const_share_e[i],
                                                                   kyber.default_thread_num),
                                             kyber.default_thread_num);
        share_U_e = kyber.pss.MultiplicationReShare(share_U_e, kyber.default_thread_num);
    }

    Vec<Vec<Pair<ZZ_p, ZZ_p>>> ret;
    ret.SetLength(2);
    ret[0] = share_U_s;
    ret[1] = share_U_e;
    return ret;
}

