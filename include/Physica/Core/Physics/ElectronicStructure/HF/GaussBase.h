/*
 * Copyright 2021-2026 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Math/Calculus/SpecialFunctions.h"
#include "Physica/Core/Scalar/Real.h"

namespace Physica {
    /**
     * Reference:
     * [1] Journal of Computational Physics 1(2), 223-244 (1966); https://doi.org/10.1016/0021-9991(66)90004-0
     * [2] Saunders V R. An Introduction to Molecular Integral Evaluation[M]. Springer Netherlands, 1975.
     */
    template<Scalar T>
    class GaussBase {
    private:
        Vector3D<T> center;
        T alpha;
        size_t l;
        size_t m;
        size_t n;
    public:
        GaussBase() = default;
        GaussBase(const Vector3D<T> center_, const T& alpha_, size_t l_, size_t m_, size_t n_);
        GaussBase(const GaussBase& base) = default;
        GaussBase(GaussBase&& base) noexcept = default;
        ~GaussBase() = default;
        /* Operators */
        GaussBase& operator=(const GaussBase& base) = default;
        GaussBase& operator=(GaussBase&& base) noexcept = default;
        /* Getters */
        [[nodiscard]] static T overlap(const GaussBase& baseA, const GaussBase& baseB);
        [[nodiscard]] static T kinetic(const GaussBase& base1, const GaussBase& base2);
        [[nodiscard]] static T nuclearAttraction(const GaussBase& base1,
                                                 const GaussBase& base2,
                                                 const Vector3D<T>& corePos);
        [[nodiscard]] static T electronRepulsion(const GaussBase& base1,
                                                 const GaussBase& base2,
                                                 const GaussBase& base3,
                                                 const GaussBase& base4);
    private:
        [[nodiscard]] T squaredNorm() const;
        [[nodiscard]] static T overlapImpl(T elemPA, T elemPB, T alpha_sum, size_t index1, size_t index2);
        [[nodiscard]] static T attractionHelper(size_t i,
                                                size_t index1,
                                                size_t index2,
                                                const T& element_pa,
                                                const T& element_pb,
                                                const T& element_cp,
                                                const T& alpha_sum);
        [[nodiscard]] static T attractionHelperG(size_t L,
                                                 size_t index1,
                                                 size_t index2,
                                                 const T& element_pa,
                                                 const T& element_pb,
                                                 const T& epsilon);
        [[nodiscard]] static T attractionHelperH(size_t i,
                                                 size_t lambda,
                                                 size_t index1,
                                                 size_t index2,
                                                 const T& element_pa,
                                                 const T& element_pb,
                                                 const T& epsilon);
        [[nodiscard]] static T repulsionHelper(size_t i,
                                               size_t index1,
                                               size_t index2,
                                               size_t index3,
                                               size_t index4,
                                               const T& element_pq,
                                               const T& element_pa,
                                               const T& element_pb,
                                               const T& element_qc,
                                               const T& element_qd,
                                               const T& epsilon1,
                                               const T& epsilon2,
                                               const T& delta);
        [[nodiscard]] static T repulsionHelperH(size_t L,
                                                size_t index1,
                                                size_t index2,
                                                const T& element1,
                                                const T& element2,
                                                const T& epsilon,
                                                bool type);
        /* Static members */
        [[nodiscard]] static T helper_f(size_t j, size_t l, size_t m, T a, T b);
        [[nodiscard]] static T helper_F(size_t v, const T& t);

        friend class Physica::Test;
    };

    template<Scalar T>
    GaussBase<T>::GaussBase(const Vector3D<T> center_, const T& alpha_, size_t l_, size_t m_, size_t n_)
            : center(center_)
            , alpha(alpha_)
            , l(l_)
            , m(m_)
            , n(n_) {
        assert(alpha.isPositive());
        /**
         * Implementation of factorial use tables when index is less than 16, it is considered 16 should be enough.
         */
        assert(l < size_t(16));
        assert(m < size_t(16));
        assert(n < size_t(16));
    }
    /**
     * Eq. (3.5) of [1]
     */
    template<Scalar T>
    T GaussBase<T>::overlap(const GaussBase& baseA, const GaussBase& baseB) {
        const T alpha_sum = baseA.alpha + baseB.alpha;
        const T inv_alpha_sum = reciprocal(alpha_sum);
        const T temp = MathConst<T>::pi * inv_alpha_sum;
        const T factor = temp * sqrt(temp);

        const T temp1 = baseA.alpha * inv_alpha_sum;
        const T temp2 = T(1) - temp1;
        const T factor2 = exp(-temp1 * baseB.alpha * (baseA.center - baseB.center).squaredNorm());

        const Vector3D<T> vecP = temp1 * baseA.center + temp2 * baseB.center;
        const Vector3D<T> vecPA = vecP - baseA.center;
        const Vector3D<T> vecPB = vecP - baseB.center;
        const T factor3_x = overlapImpl(vecPA[0], vecPB[0], alpha_sum, baseA.l, baseB.l);
        const T factor3_y = overlapImpl(vecPA[1], vecPB[1], alpha_sum, baseA.m, baseB.m);
        const T factor3_z = overlapImpl(vecPA[2], vecPB[2], alpha_sum, baseA.n, baseB.n);
        return factor * factor2 * factor3_x * factor3_y * factor3_z;
    }
    /**
     * Implement operator $-\frac{1}{2} \nabla^2$
     */
    template<Scalar T>
    T GaussBase<T>::kinetic(const GaussBase& base1, const GaussBase& base2) {
        GaussBase copy = base2;
        T result = base2.alpha * T(2 * (base2.l + base2.m + base2.n) + 3) * overlap(base1, base2);
        {
            copy.l += 2;
            const T temp1 = overlap(base1, copy);
            copy.l -= 2;
            copy.m += 2;
            const T temp2 = overlap(base1, copy);
            copy.m -= 2;
            copy.n += 2;
            const T temp3 = overlap(base1, copy);
            copy.n -= 2;
            result -= T(2) * square(base2.alpha) * (temp1 + temp2 + temp3);
        }
        if (base2.l >= 2) {
            copy.l -= 2;
            result -= T(0.5) * T(base2.l * (base2.l - 1)) * overlap(base1, copy);
            copy.l += 2;
        }
        if (base2.m >= 2) {
            copy.m -= 2;
            result -= T(0.5) * T(base2.m * (base2.m - 1)) * overlap(base1, copy);
            copy.m += 2;
        }
        if (base2.n >= 2) {
            copy.n -= 2;
            result -= T(0.5) * T(base2.n * (base2.n - 1)) * overlap(base1, copy);
        }
        return result;
    }
    /**
     * Implement operator $\frac{1}{r_c}$, where $r_c$ is vector to nuclear core.
     */
    template<Scalar T>
    T GaussBase<T>::nuclearAttraction(const GaussBase& base1,
                                      const GaussBase& base2,
                                      const Vector3D<T>& corePos) {
        const T alpha_sum = base1.alpha + base2.alpha;
        const T inv_alpha_sum = reciprocal(alpha_sum);
        const T factor = T(2) * T(M_PI) / alpha_sum;
        const T temp1 = base1.alpha * inv_alpha_sum;
        const T temp2 = T(1) - temp1;
        const T factor2 = exp(-temp1 * base2.alpha * (base1.center - base2.center).squaredNorm());

        const Vector3D<T> vector_p = temp1 * base1.center + temp2 * base2.center;
        const Vector3D<T> vector_pa = vector_p - base1.center;
        const Vector3D<T> vector_pb = vector_p - base2.center;
        const Vector3D<T> vector_cp = corePos - vector_p;
        const T temp = alpha_sum * vector_cp.squaredNorm();
        T factor3 = T(0);
        for (size_t i = 0; i <= base1.l + base2.l; ++i) {
            for (size_t j = 0; j <= base1.m + base2.m; ++j) {
                for (size_t k = 0; k <= base1.n + base2.n; ++k) {
                    factor3 += attractionHelper(i, base1.l, base2.l, vector_pa[0], vector_pb[0], vector_cp[0], alpha_sum)
                             * attractionHelper(j, base1.m, base2.m, vector_pa[1], vector_pb[1], vector_cp[1], alpha_sum)
                             * attractionHelper(k, base1.n, base2.n, vector_pa[2], vector_pb[2], vector_cp[2], alpha_sum)
                             * helper_F(i + j + k, temp);
                }
            }
        }
        return factor * factor2 * factor3;
    }

    template<Scalar T>
    T GaussBase<T>::electronRepulsion(const GaussBase& base1,
                                      const GaussBase& base2,
                                      const GaussBase& base3,
                                      const GaussBase& base4) {
        const T alpha_sum1 = base1.alpha + base3.alpha;
        const T alpha_sum2 = base2.alpha + base4.alpha;

        const T factor = T(2) * square(T(M_PI)) * sqrt(T(M_PI))
                       / (alpha_sum1 * alpha_sum2 * sqrt(alpha_sum1 + alpha_sum2));

        const T inv_alpha_sum1 = reciprocal(alpha_sum1);
        const T temp1 = base1.alpha * inv_alpha_sum1;
        const T temp2 = T(1) - temp1;
        const T factor1 = exp(-temp1 * base3.alpha * (base1.center - base3.center).squaredNorm());

        const T inv_alpha_sum2 = reciprocal(alpha_sum2);
        const T temp3 = base2.alpha * inv_alpha_sum2;
        const T temp4 = T(1) - temp3;
        const T factor2 = exp(-temp3 * base4.alpha * (base2.center - base4.center).squaredNorm());

        const Vector3D<T> vector_p = temp1 * base1.center + temp2 * base3.center;
        const Vector3D<T> vector_q = temp3 * base2.center + temp4 * base4.center;
        const Vector3D<T> vector_pq = vector_p - vector_q;
        const Vector3D<T> vector_pa = vector_p - base1.center;
        const Vector3D<T> vector_pb = vector_p - base2.center;
        const Vector3D<T> vector_qc = vector_q - base3.center;
        const Vector3D<T> vector_qd = vector_q - base4.center;

        const T epsilon1 = reciprocal(T(4) * alpha_sum1);
        const T epsilon2 = reciprocal(T(4) * alpha_sum2);
        const T delta = epsilon1 + epsilon2;
        const T temp = vector_pq.squaredNorm() / (inv_alpha_sum1 + inv_alpha_sum2);
        T result = T(0);
        for (size_t i = 0; i <= base1.l + base2.l + base3.l + base4.l; ++i) {
            for (size_t j = 0; j <= base1.m + base2.m + base3.m + base4.m; ++j) {
                for (size_t k = 0; k <= base1.n + base2.n + base3.n + base4.n; ++k) {
                    result += repulsionHelper(i, base1.l, base2.l, base3.l, base4.l, vector_pq[0], vector_pa[0], vector_pb[0], vector_qc[0], vector_qd[0], epsilon1, epsilon2, delta)
                            * repulsionHelper(j, base1.m, base2.m, base3.m, base4.m, vector_pq[1], vector_pa[1], vector_pb[1], vector_qc[1], vector_qd[1], epsilon1, epsilon2, delta)
                            * repulsionHelper(k, base1.n, base2.n, base3.n, base4.n, vector_pq[2], vector_pa[2], vector_pb[2], vector_qc[2], vector_qd[2], epsilon1, epsilon2, delta)
                            * helper_F(i + j + k, temp);
                }
            }
        }
        return factor * factor1 * factor2 * result;
    }

    template<Scalar T>
    T GaussBase<T>::squaredNorm() const {
        using Internal::doubleFactorial;
        const T temp = T(M_PI_2) / alpha;
        const T factor = temp * sqrt(temp);
        const T numerator = doubleFactorial<T>(l != 0 ? (2 * l - 1) : size_t(0))
                          * doubleFactorial<T>(m != 0 ? (2 * m - 1) : size_t(0))
                          * doubleFactorial<T>(n != 0 ? (2 * n - 1) : size_t(0));
        const T denominator = pow(T(4) * alpha, T(l + m + n));
        return factor * numerator / denominator;
    }

    template<Scalar T>
    T GaussBase<T>::overlapImpl(T elemPA, T elemPB, T alpha_sum, size_t index1, size_t index2) {
        using Internal::doubleFactorial;
        T result = T(0);
        T i_float = T(0);
        for (size_t i = 0; i <= (index1 + index2) / 2; ++i) {
            const T temp = doubleFactorial<T>(i != 0 ? (2 * i - 1) : size_t(0))
                         / pow(T(2) * alpha_sum, i_float);
            const T temp_x = helper_f(2 * i, index1, index2, elemPA, elemPB);
            result += temp_x * temp;
            i_float += T(1);
        }
        return result;
    }

    template<Scalar T>
    T GaussBase<T>::attractionHelper(size_t i,
                                     size_t index1,
                                     size_t index2,
                                     const T& element_pa,
                                     const T& element_pb,
                                     const T& element_cp,
                                     const T& alpha_sum) {
        const size_t lower = (2 * i > (index1 + index2)) ? (2 * i - index1 - index2) : size_t(0);
        const T epsilon = reciprocal(T(4) * alpha_sum);
        T result = T(0);
        for (size_t lambda = lower; lambda <= i; ++lambda)
            result += attractionHelperH(i, lambda, index1, index2, element_pa, element_pb, epsilon) * pow0(element_cp, T(lambda));
        return result;
    }
    /**
     * Implemented function $G_L$ in reference [2]
     */
    template<Scalar T>
    T GaussBase<T>::attractionHelperG(size_t L,
                                      size_t index1,
                                      size_t index2,
                                      const T& element_pa,
                                      const T& element_pb,
                                      const T& epsilon) {
        using Internal::factorial;
        T result = T(0);
        for (size_t l = 0; l <= (index1 + index2); ++l) {
            const T temp = factorial<T>(l) * helper_f(l, index1, index2, element_pa, element_pb);
            for (size_t q = 0; q <= l / 2; ++q)
                if ((l - 2 * q) == L)
                    result += temp * pow(epsilon, T(q)) / factorial<T>(q);
        }
        return result;
    }
    /**
     * Implemented function $H$ in reference [2]
     */
    template<Scalar T>
    T GaussBase<T>::attractionHelperH(size_t i,
                                      size_t lambda,
                                      size_t index1,
                                      size_t index2,
                                      const T& element_pa,
                                      const T& element_pb,
                                      const T& epsilon) {
        using Internal::factorial;
        T result = T(0);
        for (size_t L = 0; L <= (index1 + index2); ++L) {
            const T temp = attractionHelperG(L, index1, index2, element_pa, element_pb, epsilon);
            for (size_t t = 0; t <= L / 2; ++t)
                if ((L - t) == i && (L - 2 * t == lambda))
                    result += temp * pow(-epsilon, T(t)) / (factorial<T>(t) * factorial<T>(L - 2 * t));
        }
        return result;
    }

    template<Scalar T>
    T GaussBase<T>::repulsionHelper(size_t i,
                                    size_t index1,
                                    size_t index2,
                                    size_t index3,
                                    size_t index4,
                                    const T& element_pq,
                                    const T& element_pa,
                                    const T& element_pb,
                                    const T& element_qc,
                                    const T& element_qd,
                                    const T& epsilon1,
                                    const T& epsilon2,
                                    const T& delta) {
        using Internal::factorial;
        T result = T(0);
        for (size_t L1 = 0; L1 <= index1 + index2; ++L1) {
            const T factor1 = repulsionHelperH(L1, index1, index2, element_pa, element_pb, epsilon1, false);
            for (size_t L2 = 0; L2 <= index3 + index4; ++L2) {
                const T factor2 = factorial<T>(L1 + L2);
                const T factor3 = repulsionHelperH(L2, index3, index4, element_qc, element_qd, epsilon2, true);
                for (size_t t = 0; t <= (L1 + L2) / 2; ++t) {
                    if ((L1 + L2 - t) == i) {
                        result += ((t % 2 == 0) ? factor1 : -factor1) * factor2 * factor3
                                * pow0(element_pq, T(L1 + L2 - 2 * t))
                                / (factorial<T>(t) * factorial<T>(L1 + L2 - 2 * t) * pow(delta, T(L1 + L2 - t)));
                    }
                }
            }
        }
        return result;
    }
    /**
     * Implemented function $H$ in reference [2]
     */
    template<Scalar T>
    T GaussBase<T>::repulsionHelperH(size_t L,
                                     size_t index1,
                                     size_t index2,
                                     const T& element1,
                                     const T& element2,
                                     const T& epsilon,
                                     bool type) {
        using Internal::factorial;
        T result = T(0);
        const T factor = reciprocal(factorial<T>(L));
        for (size_t l = 0; l <= index1 + index2; ++l) {
            const T temp = factorial<T>(l) * helper_f(l, index1, index2, element1, element2);
            for (size_t q = 0; q <= l / 2; ++q)
                if ((l - 2 * q) == L)
                    result += (type || (l % 2 == 0) ? temp : -temp)
                            * pow(epsilon, T(l - q))
                            / factorial<T>(q);
        }
        return result * factor;
    }

    template<Scalar T>
    T GaussBase<T>::helper_f(size_t j, size_t l, size_t m, T a, T b) {
        assert(l + m >= j);
        using Internal::factorial;
        const size_t lower = j > m ? (j - m) : 0;
        const size_t upper = std::min(j, l);
        const T const_1 = factorial<T>(l) * factorial<T>(m);
        const T inv_a = a.isZero() ? T(0) : reciprocal(a);
        T result = T(0);
        T temp1 = pow0(a, T(l - lower));
        T temp2 = pow0(b, T(m + lower - j));
        for (size_t i = lower; i <= upper; ++i) {
            const T temp = const_1 / (factorial<T>(i) * factorial<T>(l - i) * factorial<T>(j - i) * factorial<T>(m + i - j));
            result += temp * temp1 * temp2;
            temp1 *= inv_a;
            temp2 *= b;
        }
        return result;
    }
    /**
     * Optimize: Note that in the implementation of gammaP function, we have calculated pow(t, v1),
     * which is unnecessary and will cancel with the pow operation here
     */
    template<Scalar T>
    T GaussBase<T>::helper_F(size_t v, const T& t) {
        const T half = T(0.5);
        const T v1 = T(v) + half;
        if (t.isZero())
            return half / v1;
        return half * pow(t, -v1) * gamma(v1) * gammaP(v1, t);
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<GaussBase<T>> {
    public:
        using ScalarType = T;
    };
}
