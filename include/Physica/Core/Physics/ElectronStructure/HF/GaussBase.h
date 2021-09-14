/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Calculus/SpetialFunctions.h"

namespace Physica::Core::Physics {
    template<class ScalarType> class GaussBase;

    namespace Internal {
        template<class T> class Traits;

        template<class T>
        class Traits<GaussBase<T>> {
        public:
            using ScalarType = T;
        };
    }
    /**
     * Reference:
     * [1] Clementi E, Davis D R. Electronic structure of large molecular systems[J]. Journal of Computational Physics, 1966, 1(2):223-244.
     * [2] Saunders V R. An Introduction to Molecular Integral Evaluation[M]. Springer Netherlands, 1975.
     */
    template<class ScalarType>
    class GaussBase {
    private:
        Vector<ScalarType, 3> center;
        ScalarType alpha;
        size_t l;
        size_t m;
        size_t n;
    public:
        GaussBase() = default;
        GaussBase(const Vector<ScalarType, 3> center_, const ScalarType& alpha_, size_t l_, size_t m_, size_t n_);
        GaussBase(const GaussBase& base) = default;
        GaussBase(GaussBase&& base) noexcept = default;
        ~GaussBase() = default;
        /* Operators */
        GaussBase& operator=(const GaussBase& base) = default;
        GaussBase& operator=(GaussBase&& base) noexcept = default;
        /* Getters */
        [[nodiscard]] static ScalarType overlap(const GaussBase& base1, const GaussBase& base2);
        [[nodiscard]] static ScalarType kinetic(const GaussBase& base1, const GaussBase& base2);
        [[nodiscard]] static ScalarType nuclearAttraction(const GaussBase& base1,
                                                          const GaussBase& base2,
                                                          const Vector<ScalarType, 3>& corePos);
        [[nodiscard]] static ScalarType electronRepulsion(const GaussBase& base1,
                                                          const GaussBase& base2,
                                                          const GaussBase& base3,
                                                          const GaussBase& base4);
    private:
        [[nodiscard]] ScalarType squaredNorm() const;
        [[nodiscard]] static ScalarType overlapImpl(const ScalarType& element_pa, const ScalarType& element_pb, const ScalarType& alpha_sum, size_t index1, size_t index2);
        [[nodiscard]] static ScalarType attractionHelper(size_t i,
                                                         size_t index1,
                                                         size_t index2,
                                                         const ScalarType& element_pa,
                                                         const ScalarType& element_pb,
                                                         const ScalarType& element_cp,
                                                         const ScalarType& alpha_sum);
        [[nodiscard]] static inline ScalarType attractionHelperG(size_t L,
                                                                 size_t index1,
                                                                 size_t index2,
                                                                 const ScalarType& element_pa,
                                                                 const ScalarType& element_pb,
                                                                 const ScalarType& epsilon);
        [[nodiscard]] static inline ScalarType attractionHelperH(size_t i,
                                                                 size_t lambda,
                                                                 size_t index1,
                                                                 size_t index2,
                                                                 const ScalarType& element_pa,
                                                                 const ScalarType& element_pb,
                                                                 const ScalarType& epsilon);
        [[nodiscard]] static ScalarType repulsionHelper(size_t i,
                                                        size_t index1,
                                                        size_t index2,
                                                        size_t index3,
                                                        size_t index4,
                                                        const ScalarType& element_pq,
                                                        const ScalarType& element_pa,
                                                        const ScalarType& element_pb,
                                                        const ScalarType& element_qc,
                                                        const ScalarType& element_qd,
                                                        const ScalarType& epsilon1,
                                                        const ScalarType& epsilon2,
                                                        const ScalarType& delta);
        [[nodiscard]] static inline ScalarType repulsionHelperH(size_t L,
                                                                size_t index1,
                                                                size_t index2,
                                                                const ScalarType& element1,
                                                                const ScalarType& element2,
                                                                const ScalarType& epsilon,
                                                                bool type);
        /* Static members */
        [[nodiscard]] static ScalarType helper_f(size_t j, size_t l, size_t m, const ScalarType& a, const ScalarType& b);
        [[nodiscard]] static ScalarType helper_F(size_t v, const ScalarType& t);

        friend class Test;
    };

    template<class ScalarType>
    GaussBase<ScalarType>::GaussBase(const Vector<ScalarType, 3> center_, const ScalarType& alpha_, size_t l_, size_t m_, size_t n_)
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

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::overlap(const GaussBase& base1, const GaussBase& base2) {
        const ScalarType alpha_sum = base1.alpha + base2.alpha;
        const ScalarType inv_alpha_sum = reciprocal(alpha_sum);
        const ScalarType temp = ScalarType(M_PI) * inv_alpha_sum;
        const ScalarType factor = temp * sqrt(temp);

        const ScalarType temp1 = base1.alpha * inv_alpha_sum;
        const ScalarType temp2 = ScalarType::One() - temp1;
        const ScalarType factor2 = exp(-temp1 * base2.alpha * (base1.center - base2.center).squaredNorm());
        
        const Vector<ScalarType, 3> vector_p = temp1 * base1.center + temp2 * base2.center;
        const Vector<ScalarType, 3> vector_pa = vector_p - base1.center;
        const Vector<ScalarType, 3> vector_pb = vector_p - base2.center;
        const ScalarType factor3_x = overlapImpl(vector_pa[0], vector_pb[0], alpha_sum, base1.l, base2.l);
        const ScalarType factor3_y = overlapImpl(vector_pa[1], vector_pb[1], alpha_sum, base1.m, base2.m);
        const ScalarType factor3_z = overlapImpl(vector_pa[2], vector_pb[2], alpha_sum, base1.n, base2.n);
        return factor * factor2 * factor3_x * factor3_y * factor3_z;
    }
    /**
     * Implement operator $-\frac{1}{2} \nabla^2$
     */
    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::kinetic(const GaussBase& base1, const GaussBase& base2) {
        GaussBase copy = base2;
        ScalarType result = base2.alpha * ScalarType(2 * (base2.l + base2.m + base2.n) + 3) * overlap(base1, base2);
        {
            copy.l += 2;
            const ScalarType temp1 = overlap(base1, copy);
            copy.l -= 2;
            copy.m += 2;
            const ScalarType temp2 = overlap(base1, copy);
            copy.m -= 2;
            copy.n += 2;
            const ScalarType temp3 = overlap(base1, copy);
            copy.n -= 2;
            result -= ScalarType::Two() * square(base2.alpha) * (temp1 + temp2 + temp3);
        }
        if (base2.l >= 2) {
            copy.l -= 2;
            result -= ScalarType(0.5) * ScalarType(base2.l * (base2.l - 1)) * overlap(base1, copy);
            copy.l += 2;
        }
        if (base2.m >= 2) {
            copy.m -= 2;
            result -= ScalarType(0.5) * ScalarType(base2.m * (base2.m - 1)) * overlap(base1, copy);
            copy.m += 2;
        }
        if (base2.n >= 2) {
            copy.n -= 2;
            result -= ScalarType(0.5) * ScalarType(base2.n * (base2.n - 1)) * overlap(base1, copy);
        }
        return result;
    }
    /**
     * Implement operator $\frac{1}{r_c}$, where $r_c$ is vector to nuclear core.
     */
    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::nuclearAttraction(const GaussBase& base1,
                                                        const GaussBase& base2,
                                                        const Vector<ScalarType, 3>& corePos) {
        const ScalarType alpha_sum = base1.alpha + base2.alpha;
        const ScalarType inv_alpha_sum = reciprocal(alpha_sum);
        const ScalarType factor = ScalarType::Two() * ScalarType(M_PI) / alpha_sum;
        const ScalarType temp1 = base1.alpha * inv_alpha_sum;
        const ScalarType temp2 = ScalarType::One() - temp1;
        const ScalarType factor2 = exp(-temp1 * base2.alpha * (base1.center - base2.center).squaredNorm());
        
        const Vector<ScalarType, 3> vector_p = temp1 * base1.center + temp2 * base2.center;
        const Vector<ScalarType, 3> vector_pa = vector_p - base1.center;
        const Vector<ScalarType, 3> vector_pb = vector_p - base2.center;
        const Vector<ScalarType, 3> vector_cp = corePos - vector_p;
        const ScalarType temp = alpha_sum * vector_cp.squaredNorm();
        ScalarType factor3 = ScalarType::Zero();
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

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::electronRepulsion(const GaussBase& base1,
                                                        const GaussBase& base2,
                                                        const GaussBase& base3,
                                                        const GaussBase& base4) {
        const ScalarType alpha_sum1 = base1.alpha + base2.alpha;
        const ScalarType alpha_sum2 = base3.alpha + base4.alpha;

        const ScalarType factor = ScalarType::Two() * square(ScalarType(M_PI)) * sqrt(ScalarType(M_PI))
                                / (alpha_sum1 * alpha_sum2 * sqrt(alpha_sum1 + alpha_sum2));

        const ScalarType inv_alpha_sum1 = reciprocal(alpha_sum1);
        const ScalarType temp1 = base1.alpha * inv_alpha_sum1;
        const ScalarType temp2 = ScalarType::One() - temp1;
        const ScalarType factor1 = exp(-temp1 * base2.alpha * (base1.center - base2.center).squaredNorm());
        
        const ScalarType inv_alpha_sum2 = reciprocal(alpha_sum2);
        const ScalarType temp3 = base3.alpha * inv_alpha_sum2;
        const ScalarType temp4 = ScalarType::One() - temp3;
        const ScalarType factor2 = exp(-temp3 * base4.alpha * (base3.center - base4.center).squaredNorm());

        const Vector<ScalarType, 3> vector_p = temp1 * base1.center + temp2 * base2.center;
        const Vector<ScalarType, 3> vector_q = temp3 * base3.center + temp4 * base4.center;
        const Vector<ScalarType, 3> vector_pq = vector_p - vector_q;
        const Vector<ScalarType, 3> vector_pa = vector_p - base1.center;
        const Vector<ScalarType, 3> vector_pb = vector_p - base2.center;
        const Vector<ScalarType, 3> vector_qc = vector_q - base3.center;
        const Vector<ScalarType, 3> vector_qd = vector_q - base4.center;

        const ScalarType epsilon1 = reciprocal(ScalarType(4) * alpha_sum1);
        const ScalarType epsilon2 = reciprocal(ScalarType(4) * alpha_sum2);
        const ScalarType delta = epsilon1 + epsilon2;
        const ScalarType temp = vector_pq.squaredNorm() / (inv_alpha_sum1 + inv_alpha_sum2);
        ScalarType result = ScalarType::Zero();
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

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::squaredNorm() const {
        using Physica::Core::Internal::doubleFactorial;
        const ScalarType temp = ScalarType(M_PI_2) / alpha;
        const ScalarType factor = temp * sqrt(temp);
        const ScalarType numerator = doubleFactorial<ScalarType>(l != 0 ? (2 * l - 1) : size_t(0))
                                   * doubleFactorial<ScalarType>(m != 0 ? (2 * m - 1) : size_t(0))
                                   * doubleFactorial<ScalarType>(n != 0 ? (2 * n - 1) : size_t(0));
        const ScalarType denominator = pow(ScalarType(4) * alpha, ScalarType(l + m + n));
        return factor * numerator / denominator;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::overlapImpl(const ScalarType& element_pa,
                                                  const ScalarType& element_pb,
                                                  const ScalarType& alpha_sum,
                                                  size_t index1,
                                                  size_t index2) {
        using Physica::Core::Internal::doubleFactorial;
        ScalarType result = ScalarType::Zero();
        ScalarType i_float = ScalarType::Zero();
        for (size_t i = 0; i <= (index1 + index2) / 2; ++i) {
            const ScalarType temp = doubleFactorial<ScalarType>(i != 0 ? (2 * i - 1) : size_t(0))
                                   / pow(ScalarType::Two() * alpha_sum, i_float);
            const ScalarType temp_x = helper_f(2 * i, index1, index2, element_pa, element_pb);
            result += temp_x * temp;
            i_float += ScalarType::One();
        }
        return result;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::attractionHelper(size_t i,
                                                       size_t index1,
                                                       size_t index2,
                                                       const ScalarType& element_pa,
                                                       const ScalarType& element_pb,
                                                       const ScalarType& element_cp,
                                                       const ScalarType& alpha_sum) {
        const size_t lower = (2 * i > (index1 + index2)) ? (2 * i - index1 - index2) : size_t(0);
        const ScalarType epsilon = reciprocal(ScalarType(4) * alpha_sum);
        ScalarType result = ScalarType::Zero();
        for (size_t lambda = lower; lambda <= i; ++lambda)
            result += attractionHelperH(i, lambda, index1, index2, element_pa, element_pb, epsilon) * pow(element_cp, ScalarType(lambda));
        return result;
    }
    /**
     * Implemented function $G_L$ in reference [2]
     */
    template<class ScalarType>
    inline ScalarType GaussBase<ScalarType>::attractionHelperG(size_t L,
                                                               size_t index1,
                                                               size_t index2,
                                                               const ScalarType& element_pa,
                                                               const ScalarType& element_pb,
                                                               const ScalarType& epsilon) {
        using Physica::Core::Internal::factorial;
        ScalarType result = ScalarType::Zero();
        for (size_t l = 0; l <= (index1 + index2); ++l) {
            const ScalarType temp = factorial<ScalarType>(l) * helper_f(l, index1, index2, element_pa, element_pb);
            for (size_t q = 0; q <= l / 2; ++q)
                if ((l - 2 * q) == L)
                    result += temp * pow(epsilon, ScalarType(q)) / factorial<ScalarType>(q);
        }
        return result;
    }
    /**
     * Implemented function $H$ in reference [2]
     */
    template<class ScalarType>
    inline ScalarType GaussBase<ScalarType>::attractionHelperH(size_t i,
                                                               size_t lambda,
                                                               size_t index1,
                                                               size_t index2,
                                                               const ScalarType& element_pa,
                                                               const ScalarType& element_pb,
                                                               const ScalarType& epsilon) {
        using Physica::Core::Internal::factorial;
        ScalarType result = ScalarType::Zero();
        for (size_t L = 0; L <= (index1 + index2); ++L) {
            const ScalarType temp = attractionHelperG(L, index1, index2, element_pa, element_pb, epsilon);
            for (size_t t = 0; t <= L / 2; ++t)
                if ((L - t) == i && (L - 2 * t == lambda))
                    result += temp * pow(-epsilon, ScalarType(t)) / (factorial<ScalarType>(t) * factorial<ScalarType>(L - 2 * t));
        }
        return result;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::repulsionHelper(size_t i,
                                                      size_t index1,
                                                      size_t index2,
                                                      size_t index3,
                                                      size_t index4,
                                                      const ScalarType& element_pq,
                                                      const ScalarType& element_pa,
                                                      const ScalarType& element_pb,
                                                      const ScalarType& element_qc,
                                                      const ScalarType& element_qd,
                                                      const ScalarType& epsilon1,
                                                      const ScalarType& epsilon2,
                                                      const ScalarType& delta) {
        using Physica::Core::Internal::factorial;
        ScalarType result = ScalarType::Zero();
        for (size_t L1 = 0; L1 <= index1 + index2; ++L1) {
            const ScalarType factor1 = repulsionHelperH(L1, index1, index2, element_pa, element_pb, epsilon1, false);
            for (size_t L2 = 0; L2 <= index3 + index4; ++L2) {
                const ScalarType factor2 = factorial<ScalarType>(L1 + L2);
                const ScalarType factor3 = repulsionHelperH(L2, index3, index4, element_qc, element_qd, epsilon2, true);
                for (size_t t = 0; t <= (L1 + L2) / 2; ++t) {
                    if ((L1 + L2 - t) == i) {
                        result += ((t % 2 == 0) ? factor1 : -factor1) * factor2 * factor3
                                * pow(element_pq, ScalarType(L1 + L2 - 2 * t))
                                / (factorial<ScalarType>(t) * factorial<ScalarType>(L1 + L2 - 2 * t) * pow(delta, ScalarType(L1 + L2 - t)));
                    }
                }
            }
        }
        return result;
    }
    /**
     * Implemented function $H$ in reference [2]
     */
    template<class ScalarType>
    inline ScalarType GaussBase<ScalarType>::repulsionHelperH(size_t L,
                                                              size_t index1,
                                                              size_t index2,
                                                              const ScalarType& element1,
                                                              const ScalarType& element2,
                                                              const ScalarType& epsilon,
                                                              bool type) {
        using Physica::Core::Internal::factorial;
        ScalarType result = ScalarType::Zero();
        const ScalarType factor = reciprocal(factorial<ScalarType>(L));
        for (size_t l = 0; l <= index1 + index2; ++l) {
            const ScalarType temp = factorial<ScalarType>(l) * helper_f(l, index1, index2, element1, element2);
            for (size_t q = 0; q <= l / 2; ++q)
                if ((l - 2 * q) == L)
                    result += (type || (l % 2 == 0) ? temp : -temp)
                            * pow(epsilon, ScalarType(l - q))
                            / factorial<ScalarType>(q);
        }
        return result * factor;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::helper_f(size_t j, size_t l, size_t m, const ScalarType& a, const ScalarType& b) {
        using Physica::Core::Internal::factorial;
        assert(l + m >= j);
        const size_t lower = j > m ? (j - m) : 0;
        const size_t upper = std::min(j, l);

        ScalarType result = ScalarType::Zero();
        ScalarType temp1 = pow(a, ScalarType(l - lower));
        ScalarType temp2 = pow(b, ScalarType(m + lower - j));
        const ScalarType const_1 = factorial<ScalarType>(l) * factorial<ScalarType>(m);
        const ScalarType inv_a = a.isZero() ? ScalarType::Zero() : reciprocal(a);
        for (size_t i = lower; i <= upper; ++i) {
            const ScalarType temp = const_1 / (factorial<ScalarType>(i)
                                              * factorial<ScalarType>(l - i)
                                              * factorial<ScalarType>(j - i)
                                              * factorial<ScalarType>(m + i - j));
            result += temp * temp1 * temp2;
            temp1 *= inv_a;
            temp2 *= b;
        }
        return result;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::helper_F(size_t v, const ScalarType& t) {
        const ScalarType half = ScalarType(0.5);
        const ScalarType v1 = ScalarType(v) + half;
        return half * pow(t, -v1) * gamma(v1) * (ScalarType::One() - gammaQ(v1, t));
    }
}
