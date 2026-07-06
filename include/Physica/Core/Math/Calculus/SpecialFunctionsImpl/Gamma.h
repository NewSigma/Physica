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

#include "Physica/Core/Math/Calculus/Differential.h"

namespace Physica {
    namespace Internal {
        /**
         * Reference:
         * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:156
         */
        template<Scalar T>
        T incompGammaImpl1(const T& a, const T& x) {
            assert(a.isPositive() && !x.isNegative());
            assert(x < a + T(1) && "[Error]: When x > a + 1, the algorithm is slow, use the other method is better");
            const T epsilon = std::numeric_limits<T>::epsilon();
            T ap = a;
            T temp = reciprocal(a);
            T sum = temp;
            do {
                ap += T(1);
                temp *= x / ap;
                sum += temp;
            } while (abs(temp) >= abs(sum) * epsilon);
            return sum * exp(T(-x + a * ln(x) - lnGamma(a)));
        }
        /**
         * Reference:
         * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:161
         */
        template<Scalar T>
        T incompGammaImpl2(const T& a, const T& x) {
            assert(a.isPositive() && !x.isNegative());
            assert(x > a + T(1) && "[Error]: When x < a + 1, the algorithm is slow, use the other method is better");
            const T epsilon = std::numeric_limits<T>::epsilon();
            const T floatMin = T(std::numeric_limits<T>::min()) / epsilon;

            T b = x + T(1) - a;
            T c = reciprocal(T(floatMin));
            T d = reciprocal(b);
            T h = d;
            T temp;
            size_t i = 1;
            do {
                T an = -T(i) * (T(i) - a); //Possible to optimize use add instead of multiply
                ++i;
                b += T(2);
                d = an * d + b;
                T copy_d(d);
                if (abs(copy_d) < floatMin)
                    d = floatMin;

                c = an / c + b;
                T copy_c(c);
                if (abs(copy_c) < floatMin)
                    c = floatMin;
                d = reciprocal(d);
                temp = c * d;
                h *= temp;
            } while (abs(temp - T(1)) >= epsilon);
            return h * exp(T(-x + a * ln(x) - lnGamma(a)));
        }
    }
    /**
     * Return the logarithm of gamma(s). s must be positive.
     *
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:156
     * [2] Lanczos, C. 1964, SIAM Journal on Numerical Analysis, ser. B, vol. 1, pp. 86-96
     */
    template<FloatPrec Prec>
    Real<Prec> lnGamma(const Real<Prec>& x) noexcept {
        assert(x.isPositive());
        using T = Real<Prec>;
        if constexpr (Prec == Double) {
            /**
             * Double version is implemented with gamma = 6 and N = 9 [1] to make full use of precision of double
             *
             * TODO: The implementation is less precise than the implementation in STL, try to find out the difference
             *
             * Compare:
             * The implementation use different number of series term for float and double, while STL provided only one implementation,
             * this version should run faster.
             */
            constexpr static std::array<float64, 9> Coeffs{228.9344030404165, -342.8104127892456, 151.3843107005646, -20.01174920149977, 0.4619036553182262, -0.0001214195995667437, -1.535239091824004E-6, 1.102873029688190E-6, -2.202670452322396E-7};
            T temp = x + T(6.5);
            temp = fma((x + T(0.5)), -ln(temp), temp);
            T ser(1.000000000000123);
            T copy(x);
            for (T coeff : Coeffs) {
                copy += T(1);
                ser += coeff / copy;
            }
            return -temp + T(0.91893853320467274178) + ln(ser / x);
        }
        else {
            static_assert(Prec == Float, "[Error]: Not implemented");
            // float version is implemented with gamma = 3 and N = 4 [1] to make full use of precision of float
            constexpr static std::array<float32, 4> Coeffs{7.6845130, -3.284681, 0.05832037, 0.0001856071};
            T temp = x + T(3.5);
            temp = fma((x + T(0.5)), -ln(temp), temp);
            T ser(0.9999998);
            T copy(x);
            for (T coeff : Coeffs) {
                copy += T(1);
                ser += coeff / copy;
            }
            return -temp + T(0.91893853320467274178) + ln(ser / x);
        }
    }

    template<FloatPrec Prec>
    Real<Prec> gamma(const Real<Prec>& s) {
        return exp(lnGamma(s));
    }

    template<FloatPrec Prec>
    Real<Prec> beta(const Real<Prec>& s1, const Real<Prec>& s2) {
        return exp(lnGamma(s1) + lnGamma(s2) - lnGamma(s1 + s2));
    }

    template<Scalar T>
    T gammaP(const T& a, const T& x) {
        assert(a.isPositive() && !x.isNegative());
        return (x < a + T(1)) ? Internal::incompGammaImpl1(a, x) : (T(1) - Internal::incompGammaImpl2(a, x));
    }

    template<Scalar T>
    T gammaQ(const T& a, const T& x) {
        assert(a.isPositive() && !x.isNegative());
        return (x < a + T(1)) ? (T(1) - Internal::incompGammaImpl1(a, x)) : Internal::incompGammaImpl2(a, x);
    }

    template<FloatPrec Prec>
    Real<Prec> bigamma(const Real<Prec>& x, const Real<Prec>& step) noexcept {
        return Differential<Real<Prec>>::ridders(lnGamma<Prec>, x, step);
    }

    template<FloatPrec Prec>
    Real<Prec> erf(const Real<Prec>& x) {
        using T = Real<Prec>;
        T x2 = square(x);
        return (x.isNegative()) ? -gammaP(T(0.5), x2) : gammaP(T(0.5), x2);
    }

    template<Scalar T>
    T erfc(const T& x) {
        T x2 = square(x);
        return (x.isNegative()) ? (T(1) + gammaP(T(0.5), x2)) : gammaQ(T(0.5), x2);
    }

    template<FloatPrec Prec>
    Real<Prec> standardNormalDistribution(const Real<Prec>& x) {
        using T = Real<Prec>;
        return (erf(x / sqrt(T(2))) + T(1)) >> 1U;
    }
}
