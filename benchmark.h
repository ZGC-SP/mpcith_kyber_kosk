//
// Created by 高谦文 on 2023/4/21.
//

#ifndef PSS_BENCHMARK_H
#define PSS_BENCHMARK_H

#include <fstream>
#include "./src/plain_LWE_prover.h"
#include "./src/frodo_prover.h"
#include "./src/kyber_prover.h"

int PackedSecretSharingMultBenchmark() {
    long N = 6500;
    long t = 250, l = 2048;
    std::ofstream out;
    out.open("mul_gate_bench.csv", std::ios::out | std::ios::trunc);
    out << "t" << "," << "l" << "," << "N" << "," << "Time Cost" << std::endl;

    std::vector<long> threshold = {100, 100, 100, 100, 100, 200, 200, 200, 200, 200};
    std::vector<long> secret_size = {128, 256, 512, 1024, 2048, 128, 256, 512, 1024, 2048};
    std::vector<long> share_size = {1118, 1737, 2990, 5488, 10483, 1028, 1429, 2229, 3828, 7374};

//    std::vector<long> threshold = {250};
//    std::vector<long> secret_size = {640};
//    std::vector<long> share_size = {2500};

    for (int i = 0; i < share_size.size(); ++i) {
        PSS::PackedSecretSharing pss(share_size[i], threshold[i], secret_size[i]);
        Vec<ZZ_p> sec_1, sec_2;
        for (int j = 0; j < secret_size[i]; ++j) {
            ZZ_p r_1 = random_ZZ_p();
            ZZ_p r_2 = random_ZZ_p();
            sec_1.append(r_1);
            sec_2.append(r_2);
        }

        Vec<Pair<ZZ_p, ZZ_p>> shares_1, shares_2;
        Vec<ZZ_p> rec_sec;

        shares_1 = pss.Share(sec_1, 32);
        shares_2 = pss.Share(sec_2, 32);

        // start timer
        auto start_e = std::chrono::high_resolution_clock::now();
        Vec<Pair<ZZ_p, ZZ_p>> mul_shares = pss.Multiplication(shares_1, shares_2, 32);
        auto end_e = std::chrono::high_resolution_clock::now();
        // end timer

        std::cout << "Multiplication Time:\t"
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end_e - start_e).count()
                  << " ms" << std::endl;

        auto start_s = std::chrono::high_resolution_clock::now();
        Vec<Pair<ZZ_p, ZZ_p>> re_shares = pss.MultiplicationReShare(mul_shares, 32, true);
        auto end_s = std::chrono::high_resolution_clock::now();

        std::cout << "Multiplication Reshared Time:\t"
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end_s - start_s).count()
                  << " ms" << std::endl;

        out << threshold[i] << ","
            << secret_size[i] << ","
            << share_size[i] << ","
            << std::chrono::duration_cast<std::chrono::milliseconds>((end_e - start_e) + (end_s - start_s)).count()
            << std::endl;
    }

    out.close();
    return 0;
}

int LinearTransformBenchmark() {
    std::ofstream out;
    out.open("LT_gate_bench.csv", std::ios::out | std::ios::trunc);
    out << "t" << "," << "l" << "," << "N" << "," << "Linear Transformation Time" << std::endl;

    std::vector<long> threshold = {100, 100, 100, 100, 100, 200, 200, 200, 200, 200};
    std::vector<long> secret_size = {128, 256, 512, 1024, 2048, 128, 256, 512, 1024, 2048};
    std::vector<long> share_size = {1118, 1737, 2990, 5488, 10483, 1028, 1429, 2229, 3828, 7374};

    for (int i = 0; i < share_size.size(); ++i) {
        PSS::PackedSecretSharing pss(share_size[i], threshold[i], secret_size[i]);
        Vec<ZZ_p> x, r;
        for (int j = 0; j < secret_size[i]; ++j) {
            ZZ_p r_1 = random_ZZ_p();
            ZZ_p r_2 = random_ZZ_p();
            x.append(r_1);
            r.append(r_2);
        }

        Mat<ZZ_p> T;
        // generate a random T for linear transformation:
        T = random_mat_ZZ_p(secret_size[i], secret_size[i]);

        Vec<Pair<ZZ_p, ZZ_p>> shares_x, shares_r;
        Vec<ZZ_p> rec_sec;

        shares_x = pss.Share(x, 32);
        shares_r = pss.Share(r, 32);

        // start timer
        auto start = std::chrono::high_resolution_clock::now();
        // add [x] + [r]
        Vec<Pair<ZZ_p, ZZ_p>> add_shares = pss.Addition(shares_x, shares_r, 32);
        // reconstruct [x + r] -> x+r
        Vec<ZZ_p> hidden_x = pss.Reconstruct(add_shares, 32);
        // linear transformation
        Vec<ZZ_p> T_x = T * hidden_x;
        // reshare the result
        Vec<Pair<ZZ_p, ZZ_p>> re_shares = pss.Share(T_x, 32); // share as [T(x + r)]
        auto end = std::chrono::high_resolution_clock::now();
        // end timer

        std::cout << "Linear Transformation Time:\t"
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                  << " ms" << std::endl;

        out << threshold[i] << ","
            << secret_size[i] << ","
            << share_size[i] << ","
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << std::endl;

        // test:
        Vec<ZZ_p> Tr = T * r;
        // share Tr
        Vec<Pair<ZZ_p, ZZ_p>> Tr_shares = pss.Share(Tr, 32);
        // check T(x + r) - Tr
        Vec<Pair<ZZ_p, ZZ_p>> check_shares = pss.Subtraction(re_shares, Tr_shares, 32); // get [Tx]
        Vec<ZZ_p> Tx = pss.Reconstruct(check_shares, 32); // get Tx
        Vec<ZZ_p> plain_Tx = T * x;
        std::cout << "Linear Transformation Correctness:\t"
                  << (Tx == plain_Tx) << std::endl;
        // test over
    }

    out.close();
    return 0;
}

int PreprocessingBenchmark() {
    std::ofstream out;
    out.open("pre-processing_bench.csv", std::ios::out | std::ios::trunc);
    out << "t" << "," << "l" << "," << "N" << "," << "Pre-processing Time" << std::endl;

    std::vector<long> threshold = {100, 100, 100, 100, 100, 200, 200, 200, 200, 200};
    std::vector<long> secret_size = {128, 256, 512, 1024, 2048, 128, 256, 512, 1024, 2048};
    std::vector<long> share_size = {1118, 1737, 2990, 5488, 10483, 1028, 1429, 2229, 3828, 7374};

    for (int i = 0; i < share_size.size(); ++i) {
        Mat<ZZ_p> A;
        // generate a random T for linear transformation:
        A = random_mat_ZZ_p(secret_size[i], secret_size[i]);

        // start timer
        auto start = std::chrono::high_resolution_clock::now();
        PlainLWE_Prover p = PlainLWE_Prover(3, 1, threshold[i], share_size[i], A, 32);
        auto end = std::chrono::high_resolution_clock::now();
        // end timer

        std::cout << "Pre-processing Time:\t"
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                  << " ms" << std::endl;

        out << threshold[i] << ","
            << secret_size[i] << ","
            << share_size[i] << ","
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << std::endl;
    }

    out.close();
    return 0;
}

int PlainLweBatchProveBenchmark() {
    // write to file
    std::ofstream out;
    out.open("plain_lwe_batch_prove.csv", std::ios::out | std::ios::trunc);
    out << "Batch Size" << "," << "Pre-processing Time" << "," << "Proving Time" << "," << "Total Time" << std::endl;

//    std::vector<long> batch = {1, 5, 10, 50, 100, 1000};
    std::vector<long> batch = {1, 5};
//    long LweSize = 640;
    long LweSize = 2048;
    std::vector<PlainLWE> instances;
    instances.resize(batch.size());
    for (int i = 0; i < instances.size(); ++i) {
        instances[i] = GeneratePlainLweInstance(LweSize, LweSize, batch[i]);
    }

    for (int i = 0; i < instances.size(); ++i) {
        std::cout << "===== Batch size: " << batch[i] << " =====" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
//        PlainLWE_Prover p = PlainLWE_Prover(3, batch[i], 250, 6500, instances[i].A, 32);
        PlainLWE_Prover p = PlainLWE_Prover(3, batch[i], 250, 6500, instances[i].A, 32);
        auto end = std::chrono::high_resolution_clock::now();

        p.InitializeConstantShares();

        auto start_p = std::chrono::high_resolution_clock::now();
        for (int j = 0; j < batch[i]; ++j) {
            p.Prove(instances[i].s[j], instances[i].e[j], instances[i].b[j]);
        }
        auto end_p = std::chrono::high_resolution_clock::now();

        std::cout << "Pre-processing Time:\t"
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                  << " ms" << std::endl;
        std::cout << "Plain LWE Prove Time:\t"
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end_p - start_p).count()
                  << " ms" << std::endl;

        std::cout << "Total Prove Time:\t"
                  << std::chrono::duration_cast<std::chrono::milliseconds>((end - start) + (end_p - start_p)).count()
                  << " ms" << std::endl;

        out << batch[i] << ","
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << ","
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_p - start_p).count() << ","
            << std::chrono::duration_cast<std::chrono::milliseconds>((end - start) + (end_p - start_p)).count()
            << std::endl;
    }

    out.close();

    return 0;
}


int KyberProverBenchmark() {
    ZZ_pPush push(ZZ(3329));

    long n = 256, m = 8;
    long k = 2;
    long N = 1454, t = 150;

    // generate random matrix to simulate NTT
    Mat<ZZ_p> T = random_mat_ZZ_p(n, n);

    // generate A:
    Mat<ZZ_pX> A;
    A.SetDims(2, 2);
    ZZ_pX a00 = random_ZZ_pX(n);
    A.put(0,0, a00);
    ZZ_pX a01 = random_ZZ_pX(n);
    A.put(0,1, a01);
    ZZ_pX a10 = random_ZZ_pX(n);
    A.put(1,0, a10);
    ZZ_pX a11 = random_ZZ_pX(n);
    A.put(1,1, a11);

    std::cout << "Start Pre-processing" << std::endl;
    auto start_pre = std::chrono::high_resolution_clock::now();
    KyberProver kyber = KyberProver(n ,N, t, k, T, 32);
    auto end_pre = std::chrono::high_resolution_clock::now();

    std::cout << "Pre-processing Time:\t"
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_pre - start_pre).count()
              << " ms" << std::endl;

    kyber.InitializeConstantShares();

    Mat<ZZ_p> S = random_mat_ZZ_p(n, 2);
    Mat<ZZ_p> E = random_mat_ZZ_p(n, 2);
    Mat<ZZ_p> B;

    std::cout << "Start Proving" << std::endl;
    auto start_prove = std::chrono::high_resolution_clock::now();
    KyberProver::Prove(kyber, A, S, E, B);
    auto end_prove = std::chrono::high_resolution_clock::now();

    std::cout << "Kyber Prove Time:\t"
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_prove - start_prove).count()
              << " ms" << std::endl;
    return 0;
}

int FrodoProverBenchmark() {
    long prime = GenPrime_long(15);
    ZZ_pPush push((ZZ(prime)));

    long n = 640, m = 8;
    long k = 9;
    long N = 2500, t = 250;

    Mat<ZZ_p> A = random_mat_ZZ_p(n, n);

    std::cout << "Start Pre-processing" << std::endl;
    auto start_pre = std::chrono::high_resolution_clock::now();
    Frodo_Prover frodo = Frodo_Prover(n, m, k, t, N, A, 32);
    auto end_pre = std::chrono::high_resolution_clock::now();

    std::cout << "Pre-processing Time:\t"
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_pre - start_pre).count()
              << " ms" << std::endl;

    frodo.InitializeConstantShares();

    Mat<ZZ_p> S = random_mat_ZZ_p(n, m);
    Mat<ZZ_p> E = random_mat_ZZ_p(n, m);
    Mat<ZZ_p> B = A * S + E;

    std::cout << "Start Proving" << std::endl;
    auto start_prove = std::chrono::high_resolution_clock::now();
    Frodo_Prover::Prove(frodo, S, E, B);
    auto end_prove = std::chrono::high_resolution_clock::now();

    std::cout << "Frodo Prove Time:\t"
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_prove - start_prove).count()
              << " ms" << std::endl;

    return 0;
}


#endif //PSS_BENCHMARK_H
