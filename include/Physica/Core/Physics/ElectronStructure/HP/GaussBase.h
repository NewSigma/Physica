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
    private:
        [[nodiscard]] ScalarType squaredNorm() const;
        [[nodiscard]] ScalarType overlapImpl(const ScalarType& element_pa, const ScalarType& element_pb, const ScalarType& alpha_sum, size_t index1, size_t index2) const;
        /* Static members */
        [[nodiscard]] static ScalarType helper_f(size_t j, size_t l, size_t m, const ScalarType& a, const ScalarType& b);
        [[nodiscard]] static ScalarType helper_F(size_t v, const ScalarType& t);
    };

    template<class ScalarType>
    GaussBase<ScalarType>::GaussBase(const Vector<ScalarType, 3> center_, const ScalarType& alpha_, size_t l_, size_t m_, size_t n_)
            : center(center_)
            , alpha(alpha_)
            , l(l_)
            , m(m_)
            , n(n_) {
        assert(alpha.isPositive());
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
        {
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
        }
        return result;
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
                                                  size_t index2) const {
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
        return half * pow(t, -half) * gamma(v1) * gammaQ(v1, t);
    }
}
