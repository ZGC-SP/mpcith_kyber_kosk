//
// Created by 高谦文 on 2023/4/6.
//

#ifndef PSS_UTILS_H
#define PSS_UTILS_H

#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/pair.h>
#include <vector>
#include <NTL/BasicThreadPool.h>

#include <openssl/sha.h>
#include <string>
#include <sstream>

using namespace NTL;

/**
 * Lagrange interpolation
 * @param points a vector(NTL::Vec) of points (in NTL::Pair<ZZ_p, ZZ_p> type)
 * @return ZZ_pX a polynomial
 */
ZZ_pX LagrangeInterpolation(const Vec<Pair<ZZ_p, ZZ_p>> &points, long thread_num = 1);

ZZ_pX LagrangeInterpolationContinuous(const Vec<Pair<ZZ_p, ZZ_p>> &points,
                                      const Vec<ZZ_p>& factorial_table,
                                      long thread_num = 1);

/**
 * Lagrange interpolation with optimization
 * @param points
 * @param thread_num
 * @return
 */
ZZ_pX NtlInterpolation(const Vec<Pair<ZZ_p, ZZ_p>> &points, long thread_num = 1);

/**
 * Convert std::vector to NTL::Vec
 * @tparam T
 * @param vec
 * @return
 */
template <typename T>
Vec<T> StdVectorToVec(const std::vector<T> &vec) {
    Vec<T> res;
    res.SetLength((long)vec.size());
    for(auto i = 0; i < vec.size(); i++) {
        res[i] = vec[i];
    }
    return res;
}

/**
 * Convert NTL::Vec to std::vector
 * @tparam T
 * @param vec
 * @return
 */
template <typename T>
std::vector<T> VecToStdVector(const Vec<T> &vec) {
    std::vector<T> res;
    res.resize((unsigned long)vec.length());
    for (int i = 0; i < vec.length(); ++i) {
        res[i] = vec[i];
    }
    return res;
}

/**
 * Generate a short vector (v[i] \in {0, 1})
 * @param len
 * @return
 */
Vec<ZZ_p> GenerateShortVector(long len);

/**
 * Generate a ranged vector (v[i] \in [-value, value])
 * @param len
 * @param value
 * @return
 */
Vec<ZZ_p> GenerateShortRangedVector(long len, const ZZ &value);

/**
 * Generate a random vector (v[i] \in Z_p)
 * @param len
 * @return
 */
Vec<ZZ_p> GenerateRandomVector(long len, const ZZ &p);

/**
 * Generate a random matrix (M[i][j] \in Z_p)
 * @param row
 * @param col
 * @return
 */
Mat<ZZ_p> GenerateRandomMatrix(long row, long col, const ZZ &p);

/**
 *
 * @param vec
 * @return
 */
ZZ_pX VectorToPolynomial(const Vec<ZZ_p> &vec);

/**
 *
 * @param poly
 * @param size
 * @return
 */
Vec<ZZ_p> PolynomialToVector(const ZZ_pX &poly, long size);

/**
 * Convert a vector to a string
 * @param vec
 * @return
 */
std::string VectorToString(const Vec<ZZ_p> &vec);

/**
 * Get the sha256 hash of a string
 * @param str
 * @return
 */
std::string SHA256(const std::string &str);

ZZ HexStringToZZ(const std::string &hex);
ZZ HexCharToZZ(char val);

// test template
template<typename Field, typename poly_Field>
poly_Field LagrangeInterpolation(const Vec<Pair<Field, Field>> &points, long thread_num = 1) {
    SetNumThreads(thread_num);

    // save the context of algebraic field
    ZZ_pContext context;
    context.save();

    // init the result polynomial f(x)
    Vec<poly_Field> F;
    F.SetLength(points.length());
    NTL_EXEC_RANGE(points.length(), first, end)
    context.restore();
    for (auto i = first; i < end; ++i) {
        poly_Field l = poly_Field();
        SetCoeff(l, 0, 0);
        F[i] = l;
    }
    NTL_EXEC_RANGE_END

    // generate the table of Mul_(i != j) (x - x[j]), the base of L_i(x)
    poly_Field T = poly_Field();
    T.SetLength(1);
    SetCoeff(T, 0, 1);
    for (const auto &point: points) {
        poly_Field t = poly_Field();
        t.SetLength(2);
        SetCoeff(t, 0, -1 * point.a);
        SetCoeff(t, 1, 1);
        T *= t;
    }
    Vec<poly_Field> subL;
    subL.SetLength(points.length());
    NTL_EXEC_RANGE(points.length(), first, end)
    context.restore();
    for (auto i = first; i < end; ++i) {
        poly_Field t = poly_Field();
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
        Field inv = Field(1);
        for (auto j = 0; j < points.length(); j++) {
            if (i != j && points[i].a != points[j].a) {
                inv *= (points[i].a - points[j].a);
            }
        }
        poly_Field l = subL[i] / inv;
        F[i] = points[i].b * l;
    }
    NTL_EXEC_RANGE_END

    poly_Field f;
    f.SetLength(long(points.length()));
    for (auto i = 0; i < points.length(); ++i) {
        f += F[i];
    }

    // std::cout << "f(x) = " << f << std::endl;
}

#endif //PSS_UTILS_H
