//
// Created by 高谦文 on 2023/4/14.
//

#include "plain_LWE_prover.h"

PlainLWE GeneratePlainLweInstance(long m, long n, long size) {
    Mat<ZZ_p> A = random_mat_ZZ_p(m, n);
    std::vector<Vec<ZZ_p>> ss, es, bs;
    for (auto i = 0; i < size; ++i) {
        Vec<ZZ_p> s = random_vec_ZZ_p(n);
        Vec<ZZ_p> e = random_vec_ZZ_p(n);
        Vec<ZZ_p> b = A * s + e;

        ss.push_back(s);
        es.push_back(e);
        bs.push_back(b);
    }
    return PlainLWE{A, ss, es, bs};
}

// Path: src/PlainLWE_Prover.cpp -> class: PlainLWE_Prover
PlainLWE_Prover::PlainLWE_Prover(long k, long v, // security parameters, batch size
                                 long t, long N, // threshold, number of parties in PSS
                                 const Mat<ZZ_p> &mat_A,
                                 long default_thread_num) {

    SetNumThreads(default_thread_num);
    ZZ_pContext context;
    context.save();


    this->A = mat_A;
    this->m = mat_A.NumRows();
    this->n = mat_A.NumCols();
    this->k = k;
    this->v = v;
    this->pss = PSS::PackedSecretSharing(N, t, this->n);

    this->default_thread_num = default_thread_num;

    // initialize the variables
    this->vec_f_shares.SetLength(this->k + this->v + 1);
    this->vec_Af_shares.SetLength(this->k + this->v + 1);

    Vec<Vec<ZZ_p>> vec_f, vec_Af;
    vec_f.SetLength(this->k + this->v + 1);
    vec_Af.SetLength(this->k + this->v + 1);

    Mat<ZZ_p> mat_f, mat_Af;
    mat_f.SetDims(this->n, this->k + this->v + 1);
    // mat_f = random_mat_ZZ_p(this->n, this->k + this->v + 1);

    mat_Af.SetDims(this->n, this->k + this->v + 1);
    // mat_Af = this->A * mat_f;

    // generate random vector r and calculate A * r.
    for (auto i = 0; i < (this->k + this->v + 1); ++i) {
        Vec<ZZ_p> f = random_vec_ZZ_p(this->n);
        vec_f[i] = f;
        Vec<ZZ_p> Af = this->A * f;
        vec_Af[i] = Af;

        // Share the random f as r[i]
        // The point (i, f(i)) is the i-th share of f
        this->vec_f_shares[i].SetLength(this->pss.GetN());
        // convert f to poly:
        ZZ_pX poly_f = VectorToPolynomial(f);
        Vec<Pair<ZZ_p, ZZ_p>> f_shares;
        f_shares.SetLength(this->pss.GetN());

        NTL_EXEC_RANGE(this->pss.GetN(), first, end)
                        context.restore();
                        for (auto j = first; j < end; ++j) {
                            Pair<ZZ_p, ZZ_p> s = Pair<ZZ_p, ZZ_p>(ZZ_p(j), eval(poly_f, ZZ_p(j)));
                            f_shares[j] = s;
                        }
        NTL_EXEC_RANGE_END
//        Vec<Pair<ZZ_p, ZZ_p>> f_shares = this->pss.Share(f, default_thread_num);
        this->vec_f_shares[i] = f_shares;

        auto end_f = std::chrono::high_resolution_clock::now();

        auto start_Af = std::chrono::high_resolution_clock::now();
        // Share the random A * f as Ar[i]
        this->vec_Af_shares[i].SetLength(this->pss.GetN());
        Vec<Pair<ZZ_p, ZZ_p>> Af_shares = this->pss.Share(Af, default_thread_num, true);
        this->vec_Af_shares[i] = Af_shares;
    }

    // make the commitment table
    this->commitments.SetDims(this->k + this->v + 1, this->pss.GetN());
    for (int i = 0; i < this->k + this->v + 1; ++i) {
        for (int j = 0; j < this->pss.GetN(); ++j) {
            std::ostringstream oss;
            oss << this->vec_f_shares[i][j].b << this->vec_Af_shares[i][j].b;
            std::string comm = SHA256(oss.str());
            this->commitments.put(i, j, comm);
        }
    }

    /**
     * Evaluate:
     *  r = F(x) = f[0] + f[1] * x + ... + f[k + v - 1] * x^(k + v - 1)
     *  Ar = A * F(x) = A * f[0] + A * f[1] * x + ... + A * f[k + v - 1] * x^(k + v - 1)
     */

    MakeMatrix(mat_f, vec_f);
    mat_f = transpose(mat_f);

    MakeMatrix(mat_Af, vec_Af);
    mat_Af = transpose(mat_Af);

    // For cut-and-choose
    ZZ_p alpha = random_ZZ_p();

    Vec<ZZ_p> vec_alpha;
    vec_alpha.SetLength(this->k + this->v + 1);
    NTL_EXEC_RANGE(this->k + this->v + 1, first, end)
                    context.restore();
                    for (auto i = first; i < end; ++i) {
                        vec_alpha[i] = power(alpha, i);
                    }
    NTL_EXEC_RANGE_END

    this->r = mat_f * vec_alpha;
    this->Ar = mat_Af * vec_alpha;

    this->r_shares = this->pss.Share(this->r, default_thread_num, true);
    this->Ar_shares = this->pss.Share(this->Ar, default_thread_num, true);

}

void PlainLWE_Prover::Prove(Vec<ZZ_p> s, Vec<ZZ_p> e, Vec<ZZ_p> b) {
    // Phase 1: prove the linear relation
    // prover shares the secret s and noise e
    // linear prove:
    auto start_linear = std::chrono::high_resolution_clock::now();
    Vec<Pair<ZZ_p, ZZ_p>> share_s = this->pss.Share(s, this->default_thread_num, true);
    Vec<Pair<ZZ_p, ZZ_p>> share_e = this->pss.Share(e, this->default_thread_num, true);

    // prover computes the [s] + [r] and reconstruct the s + r
    Vec<Pair<ZZ_p, ZZ_p>> share_s_add_r = this->pss.Addition(share_s, this->r_shares);
    Vec<ZZ_p> vec_s_plus_r = this->pss.Reconstruct(share_s_add_r, this->default_thread_num, true);

    // prover computes the A * (s + r) and shares the result: [A * (s + r)]
    Vec<ZZ_p> vec_As_add_r = this->A * vec_s_plus_r;
    Vec<Pair<ZZ_p, ZZ_p>> share_As_add_r = this->pss.Share(vec_As_add_r, this->default_thread_num, true);

    // prover computes the [A * (s + r)] - [A * r] and gets the [A * s]
    Vec<Pair<ZZ_p, ZZ_p>> share_As = this->pss.Subtraction(share_As_add_r, this->Ar_shares, this->default_thread_num);

    // prover computes the [A * s] + [e] and reconstructs the A * s + e
    Vec<Pair<ZZ_p, ZZ_p>> share_As_add_e = this->pss.Addition(share_As, share_e, this->default_thread_num);
    auto end_linear = std::chrono::high_resolution_clock::now();

    // Phase 2: prove the elements in s and e are short
    // prover shares the [s] - [1] and [e] - [1]
    auto start_short = std::chrono::high_resolution_clock::now();
    Vec<Pair<ZZ_p, ZZ_p>> share_s_sub_1 = this->pss.Subtraction(share_s, this->const_share_s, this->default_thread_num);
    Vec<Pair<ZZ_p, ZZ_p>> share_e_sub_1 = this->pss.Subtraction(share_e, this->const_share_e, this->default_thread_num);

    // prover computes the [U_s] = [s] * [s - [1]] and [U_e] = [e] * [e - [1]]
    Vec<Pair<ZZ_p, ZZ_p>> share_U_s = this->pss.Multiplication(share_s, share_s_sub_1, this->default_thread_num);
    Vec<Pair<ZZ_p, ZZ_p>> share_U_e = this->pss.Multiplication(share_e, share_e_sub_1, this->default_thread_num);
    auto end_short = std::chrono::high_resolution_clock::now();

    // show the time cost:
    auto duration_linear = std::chrono::duration_cast<std::chrono::milliseconds>(end_linear - start_linear);
    auto duration_short = std::chrono::duration_cast<std::chrono::milliseconds>(end_short - start_short);
    std::cout << "Linear Prove: " << duration_linear.count() << " ms" << std::endl;
    std::cout << "Short Prove: " << duration_short.count() << " ms" << std::endl;

}

void PlainLWE_Prover::InitializeConstantShares() {
    Vec<ZZ_p> vec_1_s, vec_1_e;

    vec_1_s.SetLength(this->n, ZZ_p(1));
    vec_1_e.SetLength(this->m, ZZ_p(1));

    this->const_share_s = this->pss.Share(vec_1_s, default_thread_num, true);
    this->const_share_e = this->pss.Share(vec_1_e, default_thread_num, true);
}

PlainLWE_Prover::~PlainLWE_Prover() = default;
