//
// Created by 高谦文 on 2023/4/21.
//

#include "frodo_prover.h"

Frodo_Prover::Frodo_Prover(long n, long m, // size of Frodo matrix / secrets
                           long k, // security parameter
                           long t, long N, // threshold, number of parties in PSS
                           const Mat<ZZ_p> &mat_A,
                           long default_thread_num) {

    this->A = mat_A;
    this->n = n; // A: n * n
    this->m = m; // S: n * m, E: n * m, B: n * m
    this->k = k;

    this->pss = PSS::PackedSecretSharing(N, t, this->n);

    this->default_thread_num = default_thread_num;

    // invoke the pre-processing of plain-LWE
    // the batch size is the cols of S
    PlainLWE_Prover pLWE = PlainLWE_Prover(this->k, this->m, t, N, this->A, this->default_thread_num);

    this->R = pLWE.GetR();
    this->AR = pLWE.GetAr();
    this->R_shares = pLWE.GetRShares();
    this->AR_shares = pLWE.GetArShares();

}

void Frodo_Prover::InitializeConstantShares() {
    // the element value of s and e is ranged from [0, +/-12]
    long range = 12;
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

void Frodo_Prover::UnitProve(const Frodo_Prover& frodo, long index, const Vec<NTL::ZZ_p> &s, const Vec<NTL::ZZ_p> &e, const Vec<NTL::ZZ_p> &b) {
    // the element value of s and e is ranged from [0, +/-12]
    long range = 12;
    auto start_linear = std::chrono::steady_clock::now();
    Vec<Pair<ZZ_p, ZZ_p>> share_s = frodo.pss.Share(s, frodo.default_thread_num, true);
    Vec<Pair<ZZ_p, ZZ_p>> share_e = frodo.pss.Share(e, frodo.default_thread_num, true);

    // Phase 1: prove the linear relationship.
    // prover computes the [s] + [r] and reconstruct the s + r
    Vec<Pair<ZZ_p, ZZ_p>> share_s_add_r = frodo.pss.Addition(share_s, frodo.R_shares, frodo.default_thread_num);
    Vec<ZZ_p> vec_s_plus_r = frodo.pss.Reconstruct(share_s_add_r, frodo.default_thread_num, true);

    // prover computes the A * (s + r) and shares the result: [A * (s + r)]
    Vec<ZZ_p> vec_As_add_r = frodo.A * vec_s_plus_r;
    Vec<Pair<ZZ_p, ZZ_p>> share_As_add_r = frodo.pss.Share(vec_As_add_r, frodo.default_thread_num, true);

    // prover computes the [A * (s + r)] - [A * r] and gets the [A * s]
    Vec<Pair<ZZ_p, ZZ_p>> share_As = frodo.pss.Subtraction(share_As_add_r, frodo.AR_shares, frodo.default_thread_num);

    // prover computes the [A * s] + [e] and reconstructs the A * s + e
    Vec<Pair<ZZ_p, ZZ_p>> share_As_add_e = frodo.pss.Addition(share_As, share_e, frodo.default_thread_num);
    auto end_linear = std::chrono::steady_clock::now();

    // Phase 2: prove the elements in s and e are short
    auto start_short = std::chrono::steady_clock::now();
    Vec<Pair<ZZ_p, ZZ_p>> share_U_s = frodo.pss.Subtraction(share_s, frodo.const_share_s[0], frodo.default_thread_num);
    for (int i = 1; i < frodo.const_share_s.length(); ++i) {
        share_U_s = frodo.pss.Multiplication(share_U_s,
                                             frodo.pss.Subtraction(share_s, frodo.const_share_s[i],
                                                                   frodo.default_thread_num),
                                             frodo.default_thread_num);
        share_U_s = frodo.pss.MultiplicationReShare(share_U_s, frodo.default_thread_num);
    }

    Vec<Pair<ZZ_p, ZZ_p>> share_U_e = frodo.pss.Subtraction(share_e, frodo.const_share_e[0], frodo.default_thread_num);
    for (int i = 1; i < frodo.const_share_e.length(); ++i) {
        share_U_e = frodo.pss.Multiplication(share_U_e,
                                             frodo.pss.Subtraction(share_e, frodo.const_share_e[i],
                                                                   frodo.default_thread_num),
                                             frodo.default_thread_num);
        share_U_e = frodo.pss.MultiplicationReShare(share_U_e, frodo.default_thread_num);
    }
    auto end_short = std::chrono::steady_clock::now();

    // show time cost:
    auto linear_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_linear - start_linear).count();
    auto short_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_short - start_short).count();
    std::cout << "linear time: " << linear_time << " ms" << std::endl;
    std::cout << "short time: " << short_time << " ms" << std::endl;

}

void Frodo_Prover::Prove(const Frodo_Prover& frodo, const Mat<NTL::ZZ_p> &mat_S, const Mat<NTL::ZZ_p> &mat_E, const Mat<NTL::ZZ_p> &mat_B) {

    Mat<ZZ_p> S_T = transpose(mat_S);
    Mat<ZZ_p> E_T = transpose(mat_E);
    Mat<ZZ_p> B_T = transpose(mat_B);

    for (int i = 0; i < mat_S.NumCols(); ++i) {
        // timer
        auto start = std::chrono::high_resolution_clock::now();
        Frodo_Prover::UnitProve(frodo, i, S_T[i], E_T[i], B_T[i]);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "unit prove time: " << elapsed.count() << std::endl;

    }

}