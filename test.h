//
// Created by 高谦文 on 2023/4/24.
//

#ifndef PSS_TEST_H
#define PSS_TEST_H

#include <iostream>
#include "./src/packed_secret_sharing.h"

int PlainTest(long thread, bool is_fast = false) {
    std::cout << "<< Plain test with " << thread << " threads >>" << std::endl;
    std::cout << "is_fast: " << (is_fast ? "true" : "false") << std::endl;

    // real parameters
    long prime = GenPrime_long(61);
    ZZ_p::init(ZZ(prime));
    long n = 2500;
    long t = 250;
    long l = 640;

//    ZZ_p::init(ZZ(17));
//    long n = 14;
//    long t = 2;
//    long l = 3;

    PSS::PackedSecretSharing pss(n, t, l);
    Vec<ZZ_p> sec;
    for (int i = 0; i < l; ++i) {
        ZZ_p r = random_ZZ_p();
        sec.append(r);
    }

    Vec<Pair<ZZ_p, ZZ_p>> shares;
    Vec<ZZ_p> rec_sec;

    std::chrono::time_point<std::chrono::high_resolution_clock> start_share, end_share;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_recon, end_recon;

    if (is_fast) {
        start_share = std::chrono::high_resolution_clock::now();
        shares = pss.Share(sec, thread, true);
        end_share = std::chrono::high_resolution_clock::now();

        start_recon = std::chrono::high_resolution_clock::now();
        rec_sec = pss.Reconstruct(shares, thread, true);
        end_recon = std::chrono::high_resolution_clock::now();
    } else {
        start_share = std::chrono::high_resolution_clock::now();
        shares = pss.Share(sec, thread);
        end_share = std::chrono::high_resolution_clock::now();

        start_recon = std::chrono::high_resolution_clock::now();
        rec_sec = pss.Reconstruct(shares, thread);
        end_recon = std::chrono::high_resolution_clock::now();
    }



    auto duration_share = std::chrono::duration_cast<std::chrono::microseconds>(end_share - start_share);
    auto duration_recon = std::chrono::duration_cast<std::chrono::microseconds>(end_recon - start_recon);

    std::cout << "==== time ====" << std::endl;
    std::cout << "share:\t" << duration_share.count() << " microseconds" << std::endl;
    std::cout << "reconstruct:\t" << duration_recon.count() << " microseconds" << std::endl;
    std::cout << "total:\t" << duration_share.count() + duration_recon.count() << " microseconds" << std::endl;



    std::cout << "==== check ====" << std::endl;
    bool is_correct = true;
    if (rec_sec.length() != sec.length()) {
        is_correct = false;
    } else {
        for (int i = 0; i < rec_sec.length(); ++i) {
            if (rec_sec[i] != sec[i]) {
                is_correct = false;
                break;
            }
        }
    }
    std::cout << "is_correct: \t" << (is_correct ? "true" : "false") << std::endl;
    return 0;
};

int TestAddEvaluation(long thread, bool is_fast = false) {
    std::cout << "<< Test add evaluation with " << thread << " threads >>" << std::endl;
    std::cout << "is_fast: " << (is_fast ? "true" : "false") << std::endl;

    long prime = GenPrime_long(61);
    ZZ_p::init(ZZ(prime));
    long n = 4500;
    long t = 30;
    long l = 2048;

    PSS::PackedSecretSharing pss(n, t, l);
    Vec<ZZ_p> sec_1, sec_2;
    for (int i = 0; i < l; ++i) {
        ZZ_p r_1 = random_ZZ_p();
        ZZ_p r_2 = random_ZZ_p();
        sec_1.append(r_1);
        sec_2.append(r_2);
    }

    Vec<Pair<ZZ_p, ZZ_p>> shares_1, shares_2;
    Vec<Pair<ZZ_p, ZZ_p>> shares_add;
    Vec<ZZ_p> rec_sec_add;


    std::chrono::time_point<std::chrono::high_resolution_clock> start_share, end_share;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_add, end_add;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_recon, end_recon;

    if (is_fast) {
        start_share = std::chrono::high_resolution_clock::now();
        // TODO: may can be parallel ?
        shares_1 = pss.FastShare(sec_1, thread);
        shares_2 = pss.FastShare(sec_2, thread);
        end_share = std::chrono::high_resolution_clock::now();

        // add share_1 and share_2
        start_add = std::chrono::high_resolution_clock::now();
        shares_add = pss.Addition(shares_1, shares_2, thread);
        end_add = std::chrono::high_resolution_clock::now();

        start_recon = std::chrono::high_resolution_clock::now();
        rec_sec_add = pss.FastReconstruct(shares_add, thread);
        end_recon = std::chrono::high_resolution_clock::now();
    } else {
        start_share = std::chrono::high_resolution_clock::now();
        // TODO: may can be parallel ?
        shares_1 = pss.Share(sec_1, thread);
        shares_2 = pss.Share(sec_2, thread);
        end_share = std::chrono::high_resolution_clock::now();

        // add share_1 and share_2
        start_add = std::chrono::high_resolution_clock::now();
        shares_add = pss.Addition(shares_1, shares_2, thread);
        end_add = std::chrono::high_resolution_clock::now();

        // reconstruct the share_1 + share_2
        start_recon = std::chrono::high_resolution_clock::now();
        rec_sec_add = pss.Reconstruct(shares_add, thread);
        end_recon = std::chrono::high_resolution_clock::now();
    }


    auto duration_share = std::chrono::duration_cast<std::chrono::microseconds>(end_share - start_share);
    auto duration_add = std::chrono::duration_cast<std::chrono::microseconds>(end_add - start_add);
    auto duration_recon = std::chrono::duration_cast<std::chrono::microseconds>(end_recon - start_recon);

    std::cout << "==== time ====" << std::endl;
    std::cout << "share:\t" << duration_share.count() << " microseconds" << std::endl;
    std::cout << "add:\t" << duration_add.count() << " microseconds" << std::endl;
    std::cout << "reconstruct:\t" << duration_recon.count() << " microseconds" << std::endl;
    std::cout << "total:\t" << duration_share.count() + duration_add.count() + duration_recon.count() << " microseconds"
              << std::endl;

    std::cout << "==== check ====" << std::endl;
    bool is_correct = true;
    if (rec_sec_add.length() != sec_1.length()) {
        is_correct = false;
    } else {
        for (int i = 0; i < rec_sec_add.length(); ++i) {
            if (rec_sec_add[i] != sec_1[i] + sec_2[i]) {
                is_correct = false;
                break;
            }
        }
    }
    std::cout << "is_correct: \t" << (is_correct ? "true" : "false") << std::endl;
    return 0;
};

// FIXME: How to reconstruct the result of multiplication?
int TestMulEvaluation(long thread) {
    std::cout << "<<Test mul evaluation with " << thread << " threads >>" << std::endl;
    long prime = GenPrime_long(61);
    ZZ_p::init(ZZ(prime));
    long n = 4500;
    long t = 30;
    long l = 2048;

    PSS::PackedSecretSharing pss(n, t, l);
    Vec<ZZ_p> sec_1, sec_2;
    for (int i = 0; i < l; ++i) {
        ZZ_p r_1 = random_ZZ_p();
        ZZ_p r_2 = random_ZZ_p();
        sec_1.append(r_1);
        sec_2.append(r_2);
    }

    auto start_share = std::chrono::high_resolution_clock::now();
    Vec<Pair<ZZ_p, ZZ_p>> shares_1 = pss.FastShare(sec_1, thread);
    Vec<Pair<ZZ_p, ZZ_p>> shares_2 = pss.FastShare(sec_2, thread);
    auto end_share = std::chrono::high_resolution_clock::now();

    // mul share_1 and share_2
    auto start_mul = std::chrono::high_resolution_clock::now();
    Vec<Pair<ZZ_p, ZZ_p>> shares_mul = pss.Multiplication(shares_1, shares_2, thread);
    auto end_mul = std::chrono::high_resolution_clock::now();

    // reconstruct the share_1 * share_2
    auto start_recon = std::chrono::high_resolution_clock::now();
    // re-parameterize the PSS
    PSS::PackedSecretSharing upPss(2 * n, 2 * t, 2 * l);
    Vec<ZZ_p> rec_sec_mul = upPss.FastReconstruct(shares_mul, thread);
    auto end_recon = std::chrono::high_resolution_clock::now();

    // end timer
    auto duration_share = std::chrono::duration_cast<std::chrono::microseconds>(end_share - start_share);
    auto duration_mul = std::chrono::duration_cast<std::chrono::microseconds>(end_mul - start_mul);
    auto duration_recon = std::chrono::duration_cast<std::chrono::microseconds>(end_recon - start_recon);
    auto duration_total = std::chrono::duration_cast<std::chrono::microseconds>(end_recon - start_share);

    std::cout << "==== time ====" << std::endl;
    std::cout << "share:\t" << duration_share.count() << " microseconds" << std::endl;
    std::cout << "mul:\t" << duration_mul.count() << " microseconds" << std::endl;
    std::cout << "reconstruct:\t" << duration_recon.count() << " microseconds" << std::endl;
    std::cout << "total:\t" << duration_total.count() << " microseconds" << std::endl;

    std::cout << "==== check ====" << std::endl;
    bool is_correct = true;
    for (int i = 0; i < sec_1.length(); ++i) {
        if (rec_sec_mul[i] != sec_1[i] * sec_2[i]) {
            is_correct = false;
            std::cout << "rec_sec_mul[" << i << "]:\t" << rec_sec_mul[i] << std::endl;
            break;
        }
    }
    std::cout << "is_correct: \t" << (is_correct ? "true" : "false") << std::endl;
    return 0;
};

#endif //PSS_TEST_H
