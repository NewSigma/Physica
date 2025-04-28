/*
 * Copyright 2021-2024 Weibo He.
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

namespace Physica {
    template<FloatPrec Prec>
    Real<Prec> legendreP(unsigned int l, const Real<Prec>& x) {
        using T = Real<Prec>;
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
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:189
     */
    template<FloatPrec Prec>
    Real<Prec> legendreP(unsigned int l, unsigned int m, const Real<Prec>& x) {
        using T = Real<Prec>;
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

    template<Scalar T>
    Complex<T> sphericalHarmomicY(unsigned int l, int m, const T& theta, const T& phi) {
        constexpr static double pi_4 = M_PI * 4;
        const unsigned int abs_m = std::abs(m);
        assert(l >= abs_m);
        const T factorial1 = Internal::factorial<T>(l - abs_m);
        const T factorial2 = Internal::factorial<T>(l + abs_m);
        const T factor = sqrt((T(2 * l + 1) * factorial1) / (T(pi_4) * factorial2));
        const T result_module = factor * legendreP(l, abs_m, cos(theta));
        const T m_phi = T(m) * phi;
        return Complex<T>(result_module * cos(m_phi), result_module * sin(m_phi));
    }

    template<Matrix M>
    HamonicRotator<M>::HamonicRotator(const M& axisRotation) : initialMat(3, 3), hamonicRotation() {
        assert(axisRotation.getRow() == 3 && axisRotation.getCol() == 3);
        initialMat(0, 0) = axisRotation(1, 1);
        initialMat(0, 1) = -axisRotation(1, 2);
        initialMat(0, 2) = axisRotation(1, 0);
        initialMat(1, 0) = -axisRotation(2, 1);
        initialMat(1, 1) = axisRotation(2, 2);
        initialMat(1, 2) = -axisRotation(2, 0);
        initialMat(2, 0) = axisRotation(0, 1);
        initialMat(2, 1) = -axisRotation(0, 2);
        initialMat(2, 2) = axisRotation(0, 0);
        hamonicRotation = initialMat;
    }

    template<Matrix M>
    bool HamonicRotator<M>::nearByMargin(double actual, double expected) {
        double diff = actual - expected;
        if (diff < 0.0)
            diff = -diff;
        // 5 bits of error in mantissa (source of '32 *')
        return diff < 32 * std::numeric_limits<double>::epsilon();
    }

    template<Matrix M>
    void HamonicRotator<M>::nextHamonicRotation() {
        using T = M::ScalarType;
        const int l = static_cast<int>((hamonicRotation.getRow() >> 1U) + 1U);
        M result(2 * l + 1, 2 * l + 1);
        for (int m = -l; m <= l; ++m) {
            for (int n = -l; n <= l; ++n) {
                T u, v, w;
                /* Get u, v, w */ {
                    T d = T(m == 0 ? 1 : 0);
                    T denom = (std::abs(n) == l ? 2.0 * l * (2.0 * l - 1) : (l + n) * (l - n));
                    const int abs_m = std::abs(m);
                    u = sqrt(T((l + m) * (l - m)) / denom);
                    v = T(0.5) * sqrt((T(1) + d) * T((l + abs_m - 1.0) * (l + abs_m)) / denom) * (T(1) - T(2) * d);
                    w = -T(0.5) * sqrt(T((l - abs_m - 1) * (l - abs_m)) / denom) * (T(1) - d);
                }
                // The functions U, V, W are only safe to call if the coefficients
                // u, v, w are not zero
                // When the coefficient is 0, these would attempt to access matrix elements that
                // are out of bounds. The list of rotations, @r, must have the @l - 1
                // previously completed band rotations. These functions are valid for l >= 2.
                if (!nearByMargin(u.toMachine(), 0.0))
                    u *= U(m, n, l);
                if (!nearByMargin(v.toMachine(), 0.0))
                    v *= V(m, n, l);
                if (!nearByMargin(w.toMachine(), 0.0))
                    w *= W(m, n, l);

                result(m + l, n + l) = T(u + v + w);
            }
        }
        hamonicRotation = std::move(result);
    }

    template<Matrix M>
    M::ScalarType HamonicRotator<M>::getCenteredElement(size_t row, size_t col) {
        const size_t matRow = hamonicRotation.getRow();
        assert(matRow % 2U == 1U);
        const size_t offset = matRow >> 1U;
        return hamonicRotation(row + offset, col + offset);
    }

    template<Matrix M>
    M::ScalarType HamonicRotator<M>::P(int i, int a, int b, int l) {
        const int i_1 = i + 1;
        if (b == l)
            return initialMat(i_1, 2) * getCenteredElement(a, l - 1) - initialMat(i_1, 0) * getCenteredElement(a, -l + 1);
        if (b == -l)
            return initialMat(i_1, 2) * getCenteredElement(a, -l + 1) + initialMat(i_1, 0) * getCenteredElement(a, l - 1);
        return initialMat(i_1, 1) * getCenteredElement(a, b);
    }

    template<Matrix M>
    M::ScalarType HamonicRotator<M>::U(int m, int n, int l) {
        return P(0, m, n, l);
    }

    template<Matrix M>
    M::ScalarType HamonicRotator<M>::V(int m, int n, int l) {
        if (m == 0)
            return P(1, 1, n, l) + P(-1, -1, n, l);
        if (m > 0) {
            const int delta = m == 1;
            return P(1, m - 1, n, l) * std::sqrt(1 + delta) - P(-1, -m + 1, n, l) * (1 - delta);
        }
        // Note there is apparent errata in [1,4,4b] dealing with this particular
        // case. [4b] writes it should be P*(1-d)+P*(1-d)^0.5
        // [1] writes it as P*(1+d)+P*(1-d)^0.5, but going through the math by hand,
        // you must have it as P*(1-d)+P*(1+d)^0.5 to form a 2^.5 term, which
        // parallels the case where m > 0.
        const int delta = m == -1;
        return P(1, m + 1, n, l) * (1 - delta) + P(-1, -m - 1, n, l) * std::sqrt(1 + delta);
    }

    template<Matrix M>
    M::ScalarType HamonicRotator<M>::W(int m, int n, int l) {
        // whenever this happens, w is also 0 so W can be anything
        if (m == 0)
            return 0.0;
        if (m > 0)
            return P(1, m + 1, n, l) + P(-1, -m - 1, n, l);
        return P(1, m - 1, n, l) - P(-1, -m + 1, n, l);
    }
}
