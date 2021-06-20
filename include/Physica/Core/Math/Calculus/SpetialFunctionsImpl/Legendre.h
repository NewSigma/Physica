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

#include <cassert>
#include "Physica/Core/MultiPrecision/ComplexScalar.h"

namespace Physica::Core {
    namespace Internal {
        template<ScalarType type, bool errorTrack>
        Scalar<type, errorTrack> factorial(unsigned int x) {
            constexpr static int size = 16;
            static const double cache[size] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000};
            if (x < size)
                return cache[x];
            else {
                double s = cache[size - 1];
                double i_d = size;
                for (unsigned int i = size; i <= x; i++) {
                    s *= i_d;
                    i_d += 1.0;
                }
                return s;
            }
        }
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> legendreP(unsigned int l, const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        assert(abs(x) <= T(1));
        T legendre_n(1);
        if (l == 0)
            return legendre_n;

        T legendre_n_1 = x;
        for (unsigned int i = 2; i <= l; ++i) {
            const T temp = (x * legendre_n_1 * T(2 * i - 1) - legendre_n * T(i - 1)) / T(i);
            legendre_n = legendre_n_1;
            legendre_n_1 = temp;
        }
        return legendre_n_1;
    }
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.189
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> legendreP(unsigned int l, unsigned int m, const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        assert(m <= l && abs(x) <= T(1));
        //Get P^m_m
        T legendre_m_n(1);
        if (m > 0) {
            const T temp = sqrt(T(T(1) - square(x)));
            T factor(1);
            for (unsigned int i = 1; i <= m; ++i) {
                legendre_m_n *= -factor * temp;
                factor += T(2);
            }
        }
        if (l == m)
            return legendre_m_n;

        T legendre_m_n_1 = x * legendre_m_n * T(2 * m + 1);
        for (unsigned int i = m + 2; i <= l; ++i) {
            const T temp = (x * legendre_m_n_1 * T(2 * i - 1) - legendre_m_n * T(i + m - 1)) / T(i - m);
            legendre_m_n = legendre_m_n_1;
            legendre_m_n_1 = temp;
        }
        return legendre_m_n_1;
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sphericalHarmomicY(unsigned int l,
                                                int m,
                                                const Scalar<type, errorTrack>& theta,
                                                const Scalar<type, errorTrack>& phi) {
        constexpr static double pi_4 = 3.141592653589793 * 4;
        using T = Scalar<type, errorTrack>;
        const unsigned int abs_m = std::abs(m);
        assert(l >= abs_m);
        const T factorial1 = Internal::factorial<type, errorTrack>(l - abs_m);
        const T factorial2 = Internal::factorial<type, errorTrack>(l + abs_m);
        const T factor = sqrt((T(2 * l + 1) * factorial1) / (T(pi_4) * factorial2));
        const T result_module = factor * legendreP(l, abs_m, cos(theta));
        const T m_phi = T(m) * phi;
        return ComplexScalar<type, errorTrack>(result_module * cos(m_phi), result_module * sin(m_phi));
    }
}
