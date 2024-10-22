//
// Created by 高谦文 on 2023/4/6.
//


#include "utils.h"


ZZ_pX LagrangeInterpolation(const Vec<Pair<ZZ_p, ZZ_p>> &points, long thread_num) {
    /**
     * Lagrange interpolation:
     *  - given points: (x[1], y[1]), (x[2], y[2]), ..., (x[n], y[n])
     *  - f(x) = Sum_(i: 1 -> n) [y[i] * L_i(x)]
     *  - L_i(x) = Mul_(i != j) [(x - x[j]) / (x[i] - x[j])]
     */

    SetNumThreads(thread_num);

    // save the context of algebraic field
    ZZ_pContext context;
    context.save();

    // init the result polynomial f(x)
    Vec<ZZ_pX> F;
    F.SetLength(points.length());
    NTL_EXEC_RANGE(points.length(), first, end)
                    context.restore();
                    for (auto i = first; i < end; ++i) {
                        ZZ_pX l = ZZ_pX();
                        SetCoeff(l, 0, 0);
                        F[i] = l;
                    }
    NTL_EXEC_RANGE_END

    // generate the table of Mul_(i != j) (x - x[j]), the base of L_i(x)
    ZZ_pX T = ZZ_pX();
    T.SetLength(1);
    SetCoeff(T, 0, 1);
    for (const auto &point: points) {
        ZZ_pX t = ZZ_pX();
        t.SetLength(2);
        SetCoeff(t, 0, -1 * point.a);
        SetCoeff(t, 1, 1);
        T *= t;
    }
    Vec<ZZ_pX> subL;
    subL.SetLength(points.length());
    NTL_EXEC_RANGE(points.length(), first, end)
                    context.restore();
                    for (auto i = first; i < end; ++i) {
                        ZZ_pX t = ZZ_pX();
                        t.SetLength(2);
                        SetCoeff(t, 0, -1 * points[i].a);
                        SetCoeff(t, 1, 1);
                        subL[i] = (T / t);
                    }
    NTL_EXEC_RANGE_END

    // interpolation
    NTL_EXEC_RANGE(points.length(), first, end)
                    context.restore();
                    for (auto i = first; i < end; i++) {
                        ZZ_p inv = ZZ_p(1);
                        for (auto j = 0; j < points.length(); j++) {
                            if (i != j && points[i].a != points[j].a) {
                                inv *= (points[i].a - points[j].a);
                            }
                        }
                        ZZ_pX l = subL[i] / inv;
                        F[i] = points[i].b * l;
                    }
    NTL_EXEC_RANGE_END

    ZZ_pX f;
    f.SetLength(long(points.length()));
    for (auto i = 0; i < points.length(); ++i) {
        f += F[i];
    }

    return f;
}

ZZ_pX LagrangeInterpolationContinuous(const Vec<Pair<ZZ_p, ZZ_p>> &points,
                                      const Vec<ZZ_p>& factorial_table ,
                                      long thread_num) {
    /**
     * Lagrange interpolation with continuous points:
     * - given points: (x[1], y[1]), (x[2], y[2]), ..., (x[n], y[n])
     * - if x[i] = i, ..., n,
     * - then: Li(x) = (x-x[1])...(x - x[i-1]) * (x - x[i+1]) * ... * (x - x[n]) / (i-1)! * (n-i)! * (-1)^(n-i+1)
     * - O(n) time complexity
     */

    SetNumThreads(thread_num);

    // save the context of algebraic field
    ZZ_pContext context;
    context.save();

    // init the result polynomial f(x)
    Vec<ZZ_pX> F;
    ZZ_pX l = ZZ_pX();
    SetCoeff(l, 0, 0);
    F.SetLength(points.length(), l);

    // generate the table of Mul_(i != j) (x - x[j]), the base of L_i(x)
    ZZ_pX T = ZZ_pX();
    T.SetLength(1);
    SetCoeff(T, 0, 1);
    for (const auto &point: points) {
        ZZ_pX t = ZZ_pX();
        t.SetLength(2);
        SetCoeff(t, 0, -1 * point.a);
        SetCoeff(t, 1, 1);
        T *= t;
    }

    // interpolation
    NTL_EXEC_RANGE(points.length(), first, end)
                    context.restore();
                    for (auto i = first; i < end; i++) {
                        // t = (x - x[0])(x - x[1])...(x - x[i-1])(x - x[i+1])...(x - x[n])
                        ZZ_pX t = ZZ_pX();
                        t.SetLength(2);
                        SetCoeff(t, 0, -1 * points[i].a);
                        SetCoeff(t, 1, 1);
                        // inv = (i-1)! * (n-i)! * (-1)^(n-i+1)
                        ZZ_p inv;
                        if ((points.length() + 1 - i) % 2 == 0 ) {
                            inv = factorial_table[points.length() - i - 1] * factorial_table[i];
                        } else {
                            inv = -1 * factorial_table[points.length() - i - 1] * factorial_table[i];
                        }
                        ZZ_pX l = (T / t) / inv;
                        F[i] = points[i].b * l;
                    }
    NTL_EXEC_RANGE_END

    ZZ_pX f;
    f.SetLength(long(points.length()));
    for (auto i = 0; i < points.length(); ++i) {
        f += F[i];
    }

    return f;
}

ZZ_pX NtlInterpolation(const Vec<Pair<ZZ_p, ZZ_p>> &points, long thread_num) {

    // parse the points
    Vec<ZZ_p> x;
    Vec<ZZ_p> y;
    x.SetLength(points.length());
    y.SetLength(points.length());
    for (int i = 0; i < points.length(); ++i) {
        x[i] = points[i].a;
        y[i] = points[i].b;
    }

    // start timer
    auto start = std::chrono::high_resolution_clock::now();
    ZZ_pX f = interpolate(x, y);
    auto end = std::chrono::high_resolution_clock::now();

//     std::cout << "f(x) = " << f << std::endl;
    std::cout << "ntl-int time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms"
              << std::endl;
    return f;
}

Vec<ZZ_p> GenerateShortVector(long len) {
    ZZ_pPush push(ZZ(2));
    return random_vec_ZZ_p(len);
}

Vec<ZZ_p> GenerateShortRangedVector(long len, const ZZ &range) {
    ZZ_pPush push(2 * range + 1);
    Vec<ZZ_p> _v =  random_vec_ZZ_p(len);
    Vec<ZZ_p> _r;
    _r.SetLength(len, to_ZZ_p(range + 1));
    return  _v - _r;
}

Vec<ZZ_p> GenerateRandomVector(long len, const ZZ &p) {
    ZZ_pPush push(p);
    return random_vec_ZZ_p(len);
}

Mat<ZZ_p> GenerateRandomMatrix(long row, long col, const ZZ &p) {
    ZZ_pPush push(p);
    return random_mat_ZZ_p(row, col);
}

ZZ_pX VectorToPolynomial(const Vec<ZZ_p> &vec) {
    ZZ_pX poly;
    poly.SetLength(vec.length());
    for (int i = 0; i < vec.length(); ++i) {
        SetCoeff(poly, i, vec[i]);
    }
    return poly;
}

Vec<ZZ_p> PolynomialToVector(const ZZ_pX &poly, long size) {
    return VectorCopy(poly, size);
}

std::string VectorToString(const Vec<ZZ_p> &vec) {
    std::string str;
    for (const auto &i: vec) {
        std::ostringstream oss;
        oss << i;
        str += oss.str();
    }
    return str;
}

std::string SHA256(const std::string &str) {
    unsigned char hash[SHA256_DIGEST_LENGTH];
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    SHA256_Update(&sha256, str.c_str(), str.size());
    SHA256_Final(hash, &sha256);
    std::string digest;
    for (unsigned char &i: hash) {
        char buf[2];
        sprintf(buf, "%02x", i);
        digest += (char *) &buf;
    }
    return digest;
}

ZZ HexStringToZZ(const std::string &hex) {
    ZZ val;
    val = to_ZZ(0);    //initialise the value to zero
    double base = 16;
    int j = 0;
    //convert the hex string to decimal string
    for (auto i = 0; i < hex.size(); i++) {
        auto _i = hex.size() - i - 1;
        val += HexCharToZZ(hex[_i]) * (to_ZZ((pow(base, j))));
        j++;
    }
    //cout << endl << "The value in decimal is " << val << endl;
    return val;
}

ZZ HexCharToZZ(char val) {
    if (val == 'A' || val == 'a') return to_ZZ(10);
    else if (val == 'B' || val == 'b') return to_ZZ(11);
    else if (val == 'C' || val == 'c') return to_ZZ(12);
    else if (val == 'D' || val == 'd') return to_ZZ(13);
    else if (val == 'E' || val == 'e') return to_ZZ(14);
    else if (val == 'F' || val == 'f') return to_ZZ(15);
    else return to_ZZ(val - '0');
}

