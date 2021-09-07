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
        /* Getters */
        [[nodiscard]] ScalarType overlap(const GaussBase& base) const;
    private:
        [[nodiscard]] ScalarType squaredNorm() const;
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
        const Vector<ScalarType, 3> vector_pa = center - vector_p;
        const Vector<ScalarType, 3> vector_pb = base.center - vector_p;
        ScalarType factor3_x = ScalarType::Zero();
        ScalarType factor3_y = ScalarType::Zero();
        ScalarType factor3_z = ScalarType::Zero();
        ScalarType i_float = ScalarType::Zero();
        for (size_t i = 0; i <= (l + base.l) / 2; ++i) {
            const ScalarType temp2 = Internal::doubleFactorial<ScalarType>(((2 * i - 1) < 2 * i) ? (2 * i - 1) : size_t(0))
                                   / pow(ScalarType::Two() * alpha_sum, i_float);
            const ScalarType temp_x = helper_f(2 * i, l, base.l, vector_pa[0], vector_pb[0]);
            factor3_x += temp_x * temp2;
            const ScalarType temp_y = helper_f(2 * i, l, base.l, vector_pa[1], vector_pb[1]);
            factor3_y += temp_y * temp2;
            const ScalarType temp_z = helper_f(2 * i, l, base.l, vector_pa[2], vector_pb[2]);
            factor3_z += temp_z * temp2;
            i_float += ScalarType::One();
        }
        return factor * factor2 * factor3_x * factor3_y * factor3_z;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::squaredNorm() const {
        const ScalarType temp = ScalarType(M_PI_2) / alpha;
        const ScalarType factor = temp * sqrt(temp);
        const ScalarType numerator = Internal::doubleFactorial<ScalarType>(((2 * l - 1) < 2 * l) ? (2 * l - 1) : size_t(0))
                                   * Internal::doubleFactorial<ScalarType>(((2 * m - 1) < 2 * m) ? (2 * m - 1) : size_t(0))
                                   * Internal::doubleFactorial<ScalarType>(((2 * n - 1) < 2 * n) ? (2 * n - 1) : size_t(0));
        const ScalarType denominator = pow(ScalarType(4) * alpha, ScalarType(l + m + n));
        return factor * numerator / denominator;
    }

    template<class ScalarType>
    ScalarType GaussBase<ScalarType>::helper_f(size_t j, size_t l, size_t m, const ScalarType& a, const ScalarType& b) {
        const size_t lower = std::max(static_cast<size_t>(0), j - m);
        const size_t upper = std::min(j, l);
        assert(l >= lower);

        ScalarType result = ScalarType::Zero();
        ScalarType temp1 = pow(a, ScalarType(l - lower));
        ScalarType temp2 = pow(b, ScalarType(m + lower - j));
        const ScalarType const_1 = lnGamma(ScalarType(l + 1)) - lnGamma(ScalarType(m + 1));
        const ScalarType inv_a = reciprocal(a);
        for (size_t i = lower; i <= upper; ++i) {
            const ScalarType temp = const_1
                                    - lnGamma(ScalarType(i + 1))
                                    - lnGamma(ScalarType(l - i + 1))
                                    + lnGamma(ScalarType(j - i + 1))
                                    - lnGamma(ScalarType(j - i - m + 1));
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
