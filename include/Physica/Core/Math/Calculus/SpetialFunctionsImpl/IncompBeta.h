/*
 * Copyright 2023 Weibo He.
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

#include "Physica/Core/Exception/BadConvergenceException.h"

namespace Physica::Core {
    namespace Internal {
        template<ScalarOption Option>
        Scalar<Option> incompBetaImpl(const Scalar<Option>& a, const Scalar<Option>& b, const Scalar<Option>& x) {
            using ScalarType = Scalar<Option>;
            constexpr int MaxIteration = 100;
            constexpr auto epsilon = std::numeric_limits<ScalarType>::epsilon();
            constexpr auto minimum = std::numeric_limits<ScalarType>::min();
            const ScalarType smallValue(minimum / epsilon);
            const ScalarType p_ab = a + b;
            const ScalarType p_a1 = a + ScalarType(1);
            const ScalarType m_a1 = a - ScalarType(1);

            ScalarType c = 1;
            ScalarType d = ScalarType(1) - p_ab * x / p_a1;
            if (abs(d) < smallValue)
                d = smallValue;
            d = reciprocal(d);
            ScalarType result = d;
            int iteration = 1;
            for (iteration = 1; iteration <= MaxIteration; ++iteration) {
                const ScalarType i = iteration;
                const ScalarType i2 = 2 * iteration;
                ScalarType aa = i * (b - i) * x / ((m_a1 + i2) * (a + i2));
                d = ScalarType(1) + aa * d;
                if (abs(d) < smallValue)
                    d = smallValue;
                c = ScalarType(1) + aa / c;
                if (abs(c) < smallValue)
                    c = smallValue;
                d = reciprocal(d);
                result *= c * d;

                aa = -(a + i) * (p_ab + i) * x / ((p_a1 + i2) * (a + i2));
                d = ScalarType(1) + aa * d;
                if (abs(d) < smallValue)
                    d = smallValue;
                c = ScalarType(1) + aa / c;
                if (abs(c) < smallValue)
                    c = smallValue;
                d = reciprocal(d);
                const ScalarType delta = c * d;
                result *= delta;
                if (abs(delta - ScalarType(1)) < ScalarType(epsilon))
                    break;
            }

            if (iteration > MaxIteration) [[unlikely]]
                throw BadConvergenceException("Cannot compute IncompBeta in the given iterations, maybe a or b is too large");
            return result;
        }
    }

    template<ScalarOption Option>
    Scalar<Option> incompBeta(const Scalar<Option>& a, const Scalar<Option>& b, const Scalar<Option>& x) {
        using ScalarType = Scalar<Option>;
        assert(x.isPositive() && x <= ScalarType(1) && "[Error]: Invalid value");
        if (x.isZero() || x == ScalarType(1)) [[unlikely]]
            return x;
        const ScalarType factor = exp(lnGamma(a + b) - lnGamma(a) - lnGamma(b) + a * ln(x) + b * ln(ScalarType(1) - x));
        const bool flag = x * (a + b + ScalarType(2)) < (a + ScalarType(1));
        if (flag)
            return factor * Internal::incompBetaImpl(a, b, x) / a;
        else
            return ScalarType(1) - factor * Internal::incompBetaImpl(b, a, ScalarType(1) - x) / b;
    }

    template<ScalarOption Option>
    inline Scalar<Option> studentT(size_t n, const Scalar<Option>& x) {
        using ScalarType = Scalar<Option>;
        const auto n1 = ScalarType(n);
        return ScalarType(1) - incompBeta(n1 * ScalarType(0.5), ScalarType(0.5), n1 / (n1 + square(x)));
    }

    template<ScalarOption Option>
    inline Scalar<Option> distributionF(const Scalar<Option>& v1, const Scalar<Option>& v2, const Scalar<Option>& x) {
        using ScalarType = Scalar<Option>;
        return ScalarType(1) - incompBeta(v1 * ScalarType(0.5), v2 * ScalarType(0.5), v2 / (v2 + v1 * x));
    }
}
