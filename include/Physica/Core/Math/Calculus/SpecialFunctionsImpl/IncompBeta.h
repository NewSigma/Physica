/*
 * Copyright 2023-2026 Weibo He.
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

namespace Physica {
    namespace Internal {
        template<FloatPrec Prec>
        Real<Prec> incompBetaImpl(Real<Prec> a, Real<Prec> b, Real<Prec> x) {
            using T = Real<Prec>;
            constexpr int MaxIteration = 100;
            constexpr auto epsilon = std::numeric_limits<T>::epsilon();
            constexpr auto minimum = std::numeric_limits<T>::min();
            const T smallValue(minimum / epsilon);
            const T p_ab = a + b;
            const T p_a1 = a + T(1);
            const T m_a1 = a - T(1);

            T c = 1;
            T d = T(1) - p_ab * x / p_a1;
            if (abs(d) < smallValue)
                d = smallValue;
            d = reciprocal(d);
            T result = d;
            int iteration = 1;
            for (iteration = 1; iteration <= MaxIteration; ++iteration) {
                const T i = iteration;
                const T i2 = 2 * iteration;
                T aa = i * (b - i) * x / ((m_a1 + i2) * (a + i2));
                d = fma(aa, d, T(1));
                if (abs(d) < smallValue)
                    d = smallValue;
                c = T(1) + aa / c;
                if (abs(c) < smallValue)
                    c = smallValue;
                d = reciprocal(d);
                result *= c * d;

                aa = -(a + i) * (p_ab + i) * x / ((p_a1 + i2) * (a + i2));
                d = fma(aa, d, T(1));
                if (abs(d) < smallValue)
                    d = smallValue;
                c = T(1) + aa / c;
                if (abs(c) < smallValue)
                    c = smallValue;
                d = reciprocal(d);
                const T delta = c * d;
                result *= delta;
                if (abs(delta - T(1)) < T(epsilon))
                    break;
            }

            if (iteration > MaxIteration) [[unlikely]]
                throw BadConvergenceException("Cannot compute IncompBeta in the given iterations, maybe a or b is too large");
            return result;
        }
    }

    template<FloatPrec Prec>
    Real<Prec> incompBeta(Real<Prec> a, Real<Prec> b, Real<Prec> x) {
        using T = Real<Prec>;
        assert(x.isPositive() && x <= T(1) && "[Error]: Invalid value");
        if (x.isZero() || x == T(1)) [[unlikely]]
            return x;
        const T factor = exp(lnGamma(a + b) - lnGamma(a) - lnGamma(b) + a * ln(x) + b * ln(T(1) - x));
        const bool flag = x * (a + b + T(2)) < (a + T(1));
        if (flag)
            return factor * Internal::incompBetaImpl(a, b, x) / a;
        else
            return T(1) - factor * Internal::incompBetaImpl(b, a, T(1) - x) / b;
    }

    template<FloatPrec Prec>
    Real<Prec> studentT(size_t n, Real<Prec> x) {
        using T = Real<Prec>;
        const auto n1 = T(n);
        return T(1) - incompBeta(n1 * T(0.5), T(0.5), n1 / (n1 + square(x)));
    }

    template<FloatPrec Prec>
    Real<Prec> distributionF(Real<Prec> v1, Real<Prec> v2, Real<Prec> x) {
        using T = Real<Prec>;
        return T(1) - incompBeta(v1 * T(0.5), v2 * T(0.5), v2 / (v2 + v1 * x));
    }
}
