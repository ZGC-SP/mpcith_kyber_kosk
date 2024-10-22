//
// Created by 高谦文 on 2023/4/6.
//

#include "packed_secret_sharing.h"
#include <NTL/BasicThreadPool.h>

namespace PSS {
    // public methods
    // default constructor
    PackedSecretSharing::PackedSecretSharing() = default;

    PackedSecretSharing::PackedSecretSharing(long N, long t, long l) {
        if (N < (2 * (l + t + 1))) {
            std::cerr << "Error: The number of parties N is less than 2 * (t + l + 1)." << std::endl;
            exit(1);
        }
        this->N = N;
        this->t = t;
        this->l = l;
        this->d = l + t + 1;

        // build a factorial table for interpolation:
        this->factorial_table.SetLength(this->N);
        factorial_table[0] = ZZ_p(1);
        factorial_table[1] = ZZ_p(1);
        for (auto i = 2; i < this->N; ++i) {
            factorial_table[i] = factorial_table[i - 1] * ZZ_p(i);
        }

    }

    // destructor
    PackedSecretSharing::~PackedSecretSharing() = default;

    /// Core implementations
    Vec<Pair<ZZ_p, ZZ_p>> PackedSecretSharing::Share(
            const Vec<ZZ_p> &secrets,
            long thread_num,
            bool isContinuous) const {
        // set the number of threads for parallel computing
        SetNumThreads(thread_num);

        // save the context of field
        ZZ_pContext context;
        context.save();

        /**
         * Remark: The number of secrets should match the parameters of PSS object.
         * - if the number of secrets is less than l, then pad 0 to the end of secrets.
         * - if the number of secrets is larger than l, then Share() only shares the first l secrets.
         */
        Vec<ZZ_p> sec = secrets;
        if (secrets.length() < l) {
            for (auto i = secrets.length(); i < l; ++i) {
                sec.append(ZZ_p(0));
            }
        }
        // result shares (n elements)
        Vec<Pair<ZZ_p, ZZ_p>> shares;
        shares.SetLength(this->N);

        // points for interpolation
        /**
         * The structure of points[]:
         * [ (0,secret[0]),     => secret points
         *   (1,secret[1]),     => secret points
         *   ...,               ...
         *   (l-1,secret[l-1]), => secret points
         *   ---------------------------------------------
         *   (l,random),        => random points to share
         *   (l+1,random),      => random points to share
         *   ...,               => ...
         *   (l+t,random),      => random points to share
         *   (l+t+1,random) ]   => random points to share
         *
         *   The number of secret points is l.
         *   The number of random points is t + 2.
         */
        Vec<Pair<ZZ_p, ZZ_p>> points;
        points.SetLength(this->d + 1); // degree + 1 = l + t + 2

        // rearrange the secret to the first l points
        NTL_EXEC_RANGE(this->l, first, end)

                        context.restore();
                        for (auto i = first; i < end; ++i) {
                            points[i] = Pair<ZZ_p, ZZ_p>(ZZ_p(i), sec[i]);
                        }

        NTL_EXEC_RANGE_END

        // generate random values for the threshold t + 1 points
        NTL_EXEC_RANGE(this->t + 2, first, end)

                        context.restore();
                        for (auto i = first; i < end; ++i) {
                            long _i = i + this->l;
                            ZZ_p r = random_ZZ_p();
                            Pair<ZZ_p, ZZ_p> random_point = Pair<ZZ_p, ZZ_p>(ZZ_p(_i), r);
                            points[_i] = random_point;
                            shares[i] = random_point;
                        }

        NTL_EXEC_RANGE_END

        ZZ_pX f;
        if (isContinuous) {
            f = LagrangeInterpolationContinuous(points, this->factorial_table, thread_num);
        } else {
            f = LagrangeInterpolation(points, thread_num);
        }

        // evaluate f(x) from x = t + 2 to x = n - 1 as rest shares
        NTL_EXEC_RANGE(this->N - (this->t + 2), first, end)
                        context.restore();
                        for (auto i = first; i < end; ++i) {
                            long _i = i + this->t + 2;
                            ZZ_p s = eval(f, ZZ_p(long(_i + this->l)));
                            shares[_i] = Pair<ZZ_p, ZZ_p>(ZZ_p(long(_i + this->l)), s);
                        }
        NTL_EXEC_RANGE_END

        return shares;
    }

    Vec<Pair<ZZ_p, ZZ_p>>
    PackedSecretSharing::FastShare(const Vec<ZZ_p> &secrets, long thread_num) const {
        auto start = std::chrono::high_resolution_clock::now();
        // set the number of threads for parallel computing
        SetNumThreads(thread_num);

        // save the context of field
        ZZ_pContext context;
        context.save();

        Vec<ZZ_p> sec = secrets;
        if (secrets.length() < l) {
            for (auto i = secrets.length(); i < l; ++i) {
                sec.append(ZZ_p(0));
            }
        }

        // result shares (n elements)
        Vec<Pair<ZZ_p, ZZ_p>> shares;
        shares.SetLength(this->N);

        /**
         * To speed up the secret-share process, we directly set the first l coefficients of
         * the polynomial f(x) as the secrets. The rest t + 2 coefficients are random values.
         */
        Vec<ZZ_p> polynomial_coeffs;
        polynomial_coeffs.SetLength(this->d + 1); // degree + 1 = l + t + 2
//        polynomial_coeffs.insert(polynomial_coeffs.begin(), sec.begin(), sec.begin() + this->l);
        for (int i = 0; i < this->l; ++i) {
            polynomial_coeffs[i] = sec[i];
        }

        // generate random values for the threshold t + 2 points
        NTL_EXEC_RANGE(this->t + 2, first, end)

                        context.restore();
                        for (auto i = first; i < end; ++i) {
                            long _i = i + this->l;
                            ZZ_p r = random_ZZ_p();
                            polynomial_coeffs[_i] = r;
                        }

        NTL_EXEC_RANGE_END

        // Set the polynomial f(x) = a0 + a1 * x + a2 * x^2 + ... + ad * x^d, where ai is the i-th coefficient.
        ZZ_pX f;
        f.SetLength(this->d + 1);
        for (auto i = 0; i < this->d + 1; ++i) {
            SetCoeff(f, i, polynomial_coeffs[i]);
        }

        // evaluate f(x) from x = 1 to x = n - 1 as shares
        NTL_EXEC_RANGE(this->N, first, end)

                        context.restore();
                        for (auto i = first; i < end; ++i) {
                            ZZ_p s = eval(f, ZZ_p(long(i)));
                            shares[i] = Pair<ZZ_p, ZZ_p>(ZZ_p(long(i)), s);
                        }

        NTL_EXEC_RANGE_END

        auto end = std::chrono::high_resolution_clock::now();
        // std::cout << "fast share: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
        return shares;
    }

    Vec<ZZ_p> PackedSecretSharing::Reconstruct(
            const Vec<Pair<ZZ_p, ZZ_p>> &shares,
            long thread_num,
            bool isContinuous) const {
        // if the number of shares is less than l + t + 2 (degree + 1), reconstruct fails.
        if (shares.length() < this->d + 1) {
            // throw an exception for reconstruction failure
            std::cerr << "Error: The number of shares is less than l + t + 2." << std::endl;
            exit(1);
        }

        // set the number of threads for parallel computing
        SetNumThreads(thread_num);

        // save the context of field
        ZZ_pContext context;
        context.save();

        // reconstruct the polynomial f(x) using Lagrange interpolation
        Vec<Pair<ZZ_p, ZZ_p>> used_shares;
        used_shares.SetLength(this->d + 1);
        for (int i = 0; i < this->d + 1; ++i) {
            used_shares[i] = shares[i];
        }

        ZZ_pX f;
        if (isContinuous) {
            f = LagrangeInterpolationContinuous(used_shares, this->factorial_table, thread_num);
        } else {
            f = LagrangeInterpolation(used_shares, thread_num);
        }

        // evaluate f(x) from x = 0 to x = l - 1 as recovered secrets
        Vec<ZZ_p> recovered_secrets;
        recovered_secrets.SetLength(this->l);
        NTL_EXEC_RANGE(this->l, first, end)

                        context.restore();
                        for (auto i = first; i < end; ++i) {
                            recovered_secrets[i] = eval(f, ZZ_p(i));
                        }

        NTL_EXEC_RANGE_END

        return recovered_secrets;
    }

    Vec<ZZ_p> PackedSecretSharing::FastReconstruct(const Vec<Pair<ZZ_p, ZZ_p>> &shares,
                                                   long thread_num) const {
        auto start = std::chrono::high_resolution_clock::now();
        // if the number of shares is less than l + t + 2 (degree + 1), reconstruct fails.
        if (shares.length() < this->d + 1) {
            // throw an exception for reconstruction failure
            std::cerr << "Error: The number of shares is less than l + t + 2." << std::endl;
            exit(1);
        }
        // set the number of threads for parallel computing
        SetNumThreads(thread_num);

        // save the context of field
        ZZ_pContext context;
        context.save();

        // reconstruct the polynomial f(x) using Lagrange interpolation
        Vec<Pair<ZZ_p, ZZ_p>> used_shares;
        used_shares.SetLength(this->d + 1);
        for (int i = 0; i < this->d + 1; ++i) {
            used_shares[i] = shares[i];
        }

        // test NTL interpolation
        // ZZ_pX f = NtlInterpolation(used_shares, thread_num);

        ZZ_pX f = LagrangeInterpolation(used_shares, thread_num);

        // Here we just need to get the first l coefficients of f(x) as the recovered secrets
        Vec<ZZ_p> recovered_secrets;
        recovered_secrets.SetLength(this->l);
        NTL_EXEC_RANGE(this->l, first, end)

                        context.restore();
                        for (auto i = first; i < end; ++i) {
                            recovered_secrets[i] = coeff(f, i);
                        }

        NTL_EXEC_RANGE_END

        auto end = std::chrono::high_resolution_clock::now();
        // std::cout << "fast recon: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
        return recovered_secrets;
    }

    Vec<Pair<ZZ_p, ZZ_p>> PackedSecretSharing::Addition(const Vec<Pair<ZZ_p, ZZ_p>> &x,
                                                        const Vec<Pair<ZZ_p, ZZ_p>> &y,
                                                        long thread_num) const {
        auto start = std::chrono::high_resolution_clock::now();
        // Check if the length of share-1 and share-2 are equal and if the length is larger than threshold
        if (x.length() != y.length()) {
            // throw an exception
            std::cerr << "Error: The length of share-1 and share-2 are not equal." << std::endl;
            exit(1);
        }
        if (x.length() <= this->t) {
            // throw an exception
            std::cerr << "Error: The length of share-1 and share-2 are less than threshold t." << std::endl;
            exit(1);
        }

        // set the number of threads for parallel computing
        SetNumThreads(thread_num);
        // save the context of field
        ZZ_pContext context;
        context.save();

        Vec<Pair<ZZ_p, ZZ_p>> z;
        z.SetLength(x.length());

        NTL_EXEC_RANGE(x.length(), first, end)

                        context.restore();
                        for (auto i = first; i < end; ++i) {
                            if (x[i].a != y[i].a) {
                                // throw an exception
                                std::cerr << "Error: The x-coordinate of share-1 and share-2 are not equal at [" << i
                                          << "]."
                                          << std::endl;
                                exit(1);
                            }
                            z[i] = Pair<ZZ_p, ZZ_p>(x[i].a, x[i].b + y[i].b);
                        }

        NTL_EXEC_RANGE_END
        auto end = std::chrono::high_resolution_clock::now();
        // std::cout << "add: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
        return z;
    }

    Vec<Pair<ZZ_p, ZZ_p>> PackedSecretSharing::Multiplication(const Vec<Pair<ZZ_p, ZZ_p>> &x,
                                                              const Vec<Pair<ZZ_p, ZZ_p>> &y,
                                                              long thread_num) const {
        auto start = std::chrono::high_resolution_clock::now();
        // Check if the length of share-1 and share-2 are equal and if the length is larger than threshold
        if (x.length() != y.length()) {
            // throw an exception
            std::cerr << "Error: The length of share-1 and share-2 are not equal." << std::endl;
            exit(1);
        }
        if (x.length() <= this->t) {
            // throw an exception
            std::cerr << "Error: The length of share-1 and share-2 are less than threshold t." << std::endl;
            exit(1);
        }

        // set the number of threads for parallel computing
        SetNumThreads(thread_num);
        // save the context of field
        ZZ_pContext context;
        context.save();

        Vec<Pair<ZZ_p, ZZ_p>> z;
        z.SetLength(x.length());

        NTL_EXEC_RANGE(x.length(), first, end)

                        context.restore();
                        for (auto i = first; i < end; ++i) {
                            if (x[i].a != y[i].a) {
                                // throw an exception
                                std::cerr << "Error: The x-coordinate of share-1 and share-2 are not equal at [" << i
                                          << "]."
                                          << std::endl;
                                exit(1);
                            }
                            z[i] = Pair<ZZ_p, ZZ_p>(x[i].a, x[i].b * y[i].b);
                        }

        NTL_EXEC_RANGE_END
        auto end = std::chrono::high_resolution_clock::now();
        // std::cout << "mul: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
        return z;
    }

    Vec<Pair<ZZ_p, ZZ_p>>
    PackedSecretSharing::Subtraction(const Vec<Pair<ZZ_p, ZZ_p>> &x, const Vec<Pair<ZZ_p, ZZ_p>> &y,
                                     long thread_num) const {
        auto start = std::chrono::high_resolution_clock::now();
        // Check if the length of share-1 and share-2 are equal and if the length is larger than threshold
        if (x.length() != y.length()) {
            // throw an exception
            std::cerr << "Error: The length of share-1 and share-2 are not equal." << std::endl;
            exit(1);
        }
        if (x.length() <= this->t) {
            // throw an exception
            std::cerr << "Error: The length of share-1 and share-2 are less than threshold t." << std::endl;
            exit(1);
        }

        // set the number of threads for parallel computing
        SetNumThreads(thread_num);
        // save the context of field
        ZZ_pContext context;
        context.save();

        Vec<Pair<ZZ_p, ZZ_p>> z;
        z.SetLength(x.length());

        NTL_EXEC_RANGE(x.length(), first, end)

                        context.restore();
                        for (auto i = first; i < end; ++i) {
                            if (x[i].a != y[i].a) {
                                // throw an exception
                                std::cerr << "Error: The x-coordinate of share-1 and share-2 are not equal at [" << i
                                          << "]."
                                          << std::endl;
                                exit(1);
                            }
                            z[i] = Pair<ZZ_p, ZZ_p>(x[i].a, x[i].b - y[i].b);
                        }

        NTL_EXEC_RANGE_END
        auto end = std::chrono::high_resolution_clock::now();
        // std::cout << "sub: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << std::endl;
        return z;
    }

    Vec<Pair<ZZ_p, ZZ_p>>
    PackedSecretSharing::MultiplicationReShare(const Vec<Pair<ZZ_p, ZZ_p>> &mul_shares, long thread_num, bool isContinuous) const {
        // share size: N
        // the degree of reconstructed-polynomial: 2 * d = 2 * (t + l + 1)
        // required interpolation points: 2 * d + 1 = 2 * (t + l + 1) + 1
        // if the number of shares is less than l + t + 2 (degree + 1), reconstruct fails.
        if (mul_shares.length() < (2 * this->d) + 1) {
            // throw an exception for reconstruction failure
            std::cerr << "Error: The number of shares is less than 2 * (t + l + 1) + 1." << std::endl;
            exit(1);
        }

        // set the number of threads for parallel computing
        SetNumThreads(thread_num);

        // save the context of field
        ZZ_pContext context;
        context.save();

        // reconstruct the polynomial f(x) using Lagrange interpolation
        Vec<Pair<ZZ_p, ZZ_p>> used_shares;
        used_shares.SetLength((2 * this->d) + 1);
        for (int i = 0; i < (2 * this->d) + 1; ++i) {
            used_shares[i] = mul_shares[i];
        }

        ZZ_pX f;
        if (isContinuous) {
            f = LagrangeInterpolationContinuous(used_shares, this->factorial_table, thread_num);
        } else {
            f = LagrangeInterpolation(used_shares, thread_num);
        }

        // evaluate f(x) from x = 0 to x = l - 1 as recovered secrets
        Vec<ZZ_p> recovered_secrets;
        recovered_secrets.SetLength(this->l);
        NTL_EXEC_RANGE(this->l, first, end)
                        context.restore();
                        for (auto i = first; i < end; ++i) {
                            recovered_secrets[i] = eval(f, ZZ_p(i));
                        }
        NTL_EXEC_RANGE_END

//        // debug
//        std::cout << "reconstructed shares: " << recovered_secrets << std::endl;

        // share above secrets ([A * B])
        Vec<Pair<ZZ_p, ZZ_p>> shares = this->Share(recovered_secrets, thread_num, true);
        return shares;
    }

} // PSS