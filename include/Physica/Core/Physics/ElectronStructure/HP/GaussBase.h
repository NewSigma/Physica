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
    /**
     * Reference:
     * [1] Clementi E, Davis D R. Electronic structure of large molecular systems[J]. Journal of Computational Physics, 1966, 1(2):223-244.
     * [2] Saunders V R. An Introduction to Molecular Integral Evaluation[M]. Springer Netherlands, 1975.
     */
    template<class ScalarType>
    class GaussBase {
        Vector<ScalarType, 3> center;
        ScalarType alpha;
        size_t l;
        size_t m;
        size_t n;
    public:
        GaussBase(const Vector<ScalarType, 3> center_, const ScalarType& alpha_, size_t l_, size_t m_, size_t n_);
        GaussBase(const GaussBase& base) = default;
        GaussBase(GaussBase&& base) noexcept = default;
        ~GaussBase() = default;
        /* Operators */
        GaussBase& operator=(const GaussBase& base) = default;
        GaussBase& operator=(GaussBase&& base) noexcept = default;
        /* Getters */
        [[nodiscard]] ScalarType overlap(const GaussBase& base) const;
        [[nodiscard]] ScalarType kinetic(const GaussBase& base) const;
        [[nodiscard]] ScalarType nuclearAttraction(const GaussBase& base, const Vector<ScalarType, 3>& corePos) const;
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
    ScalarType GaussBase<ScalarType>::overlap(const GaussBase& base) const {
        const ScalarType alpha_sum = alpha + base.alpha;
        const ScalarType inv_alpha_sum = reciprocal(alpha_sum);
        const ScalarType temp = ScalarType(M_PI) * inv_alpha_sum;
        const ScalarType factor = temp * sqrt(temp);

        const ScalarType temp1 = alpha * inv_alpha_sum;
        const ScalarType temp2 = ScalarType::One() - temp1;
        const ScalarType factor2 = exp(-temp1 * base.alpha * (center - base.center).squaredNorm());
        
        const Vector<ScalarType, 3> vector_p = temp1 * center + temp2 * base.center;
        const Vector<ScalarType, 3> vector_pa = vector_p - center;
        const Vector<ScalarType, 3> vector_pb = vector_p - base.center;
        const ScalarType factor3_x = overlapImpl(vector_pa[0], vector_pb[0], alpha_sum, l, base.l);
        const ScalarType factor3_y = overlapImpl(vector_pa[1], vector_pb[1], alpha_sum, m, base.m);
        const ScalarType factor3_z = overlapImpl(vector_pa[2], vector_pb[2], alpha_sum, n, base.n);
        return factor * factor2 * factor3_x * factor3_y * factor3_z;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::kinetic(const GaussBase& base) const {
        GaussBase copy = base;
        ScalarType result = base.alpha * ScalarType(2 * (base.l + base.m + base.n) + 3) * overlap(base);
        {
            copy.l += 2;
            const ScalarType temp1 = overlap(copy);
            copy.l -= 2;
            copy.m += 2;
            const ScalarType temp2 = overlap(copy);
            copy.m -= 2;
            copy.n += 2;
            const ScalarType temp3 = overlap(copy);
            copy.n -= 2;
            result -= ScalarType::Two() * square(base.alpha) * (temp1 + temp2 + temp3);
        }
        if (base.l >= 2) {
            copy.l -= 2;
            result -= ScalarType(0.5) * ScalarType(base.l * (base.l - 1)) * overlap(copy);
            copy.l += 2;
        }
        if (base.m >= 2) {
            copy.m -= 2;
            result -= ScalarType(0.5) * ScalarType(base.m * (base.m - 1)) * overlap(copy);
            copy.m += 2;
        }
        if (base.n >= 2) {
            copy.n -= 2;
            result -= ScalarType(0.5) * ScalarType(base.n * (base.n - 1)) * overlap(copy);
        }
        return result;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::nuclearAttraction(const GaussBase& base, const Vector<ScalarType, 3>& corePos) const {
        const ScalarType alpha_sum = alpha + base.alpha;
        const ScalarType inv_alpha_sum = reciprocal(alpha_sum);
        const ScalarType factor = ScalarType::Two() * ScalarType(M_PI) / alpha_sum;
        const ScalarType temp1 = alpha * inv_alpha_sum;
        const ScalarType temp2 = ScalarType::One() - temp1;
        const ScalarType factor2 = exp(-temp1 * base.alpha * (center - base.center).squaredNorm());
        
        const Vector<ScalarType, 3> vector_p = temp1 * center + temp2 * base.center;
        const Vector<ScalarType, 3> vector_pa = vector_p - center;
        const Vector<ScalarType, 3> vector_pb = vector_p - base.center;
        const Vector<ScalarType, 3> vector_cp = corePos - vector_p;
        const ScalarType temp = alpha_sum * vector_cp.squaredNorm();
        ScalarType factor3 = ScalarType::Zero();
        for (size_t i = 0; i <= l + base.l; ++i) {
            for (size_t j = 0; j <= m + base.m; ++j) {
                for (size_t k = 0; k <= n + base.n; ++k) {
                    factor3 += attractionHelper(i, l, base.l, vector_pa[0], vector_pb[0], vector_cp[0], alpha_sum)
                             * attractionHelper(j, m, base.m, vector_pa[1], vector_pb[1], vector_cp[1], alpha_sum)
                             * attractionHelper(k, n, base.n, vector_pa[2], vector_pb[2], vector_cp[2], alpha_sum)
                             * helper_F(i + j + k, temp);
                }
            }
        }
        return factor * factor2 * factor3;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::squaredNorm() const {
        const ScalarType temp = ScalarType(M_PI_2) / alpha;
        const ScalarType factor = temp * sqrt(temp);
        const ScalarType numerator = Internal::doubleFactorial<ScalarType>(l != 0 ? (2 * l - 1) : size_t(0))
                                   * Internal::doubleFactorial<ScalarType>(m != 0 ? (2 * m - 1) : size_t(0))
                                   * Internal::doubleFactorial<ScalarType>(n != 0 ? (2 * n - 1) : size_t(0));
        const ScalarType denominator = pow(ScalarType(4) * alpha, ScalarType(l + m + n));
        return factor * numerator / denominator;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::overlapImpl(const ScalarType& element_pa,
                                                  const ScalarType& element_pb,
                                                  const ScalarType& alpha_sum,
                                                  size_t index1,
                                                  size_t index2) {
        ScalarType result = ScalarType::Zero();
        ScalarType i_float = ScalarType::Zero();
        for (size_t i = 0; i <= (index1 + index2) / 2; ++i) {
            const ScalarType temp = Internal::doubleFactorial<ScalarType>(i != 0 ? (2 * i - 1) : size_t(0))
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
        ScalarType result = ScalarType::Zero();
        for (size_t l = 0; l <= (index1 + index2); ++l) {
            const ScalarType temp = Internal::factorial<ScalarType>(l) * helper_f(l, index1, index2, element_pa, element_pb);
            for (size_t q = 0; q <= l / 2; ++q)
                if ((l - 2 * q) == L)
                    result += temp * pow(epsilon, ScalarType(q)) / Internal::factorial<ScalarType>(q);
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
        ScalarType result = ScalarType::Zero();
        for (size_t L = 0; L <= (index1 + index2); ++L) {
            const ScalarType temp = attractionHelperG(L, index1, index2, element_pa, element_pb, epsilon);
            for (size_t t = 0; t <= L / 2; ++t)
                if ((L - t) == i && (L - 2 * t == lambda))
                    result += temp * pow(-epsilon, ScalarType(t)) / (Internal::factorial<ScalarType>(t) * Internal::factorial<ScalarType>(L - 2 * t));
        }
        return result;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::helper_f(size_t j, size_t l, size_t m, const ScalarType& a, const ScalarType& b) {
        assert(l + m >= j);
        const size_t lower = j > m ? (j - m) : 0;
        const size_t upper = std::min(j, l);

        ScalarType result = ScalarType::Zero();
        ScalarType temp1 = pow(a, ScalarType(l - lower));
        ScalarType temp2 = pow(b, ScalarType(m + lower - j));
        const ScalarType const_1 = lnGamma(ScalarType(l + 1)) + lnGamma(ScalarType(m + 1));
        const ScalarType inv_a = reciprocal(a);
        for (size_t i = lower; i <= upper; ++i) {
            const ScalarType temp = const_1
                                    - lnGamma(ScalarType(i + 1))
                                    - lnGamma(ScalarType(l - i + 1))
                                    - lnGamma(ScalarType(j - i + 1))
                                    - lnGamma(ScalarType(m - j - i + 1));
            result += exp(temp) * temp1 * temp2;
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
