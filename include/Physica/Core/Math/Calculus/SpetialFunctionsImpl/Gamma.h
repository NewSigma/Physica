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

#include "Physica/Core/Math/Calculus/Differential.h"

namespace Physica::Core {
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
    template<Scalar T>
    T lnGamma(const T& s) {
        assert(s.isPositive());
        if constexpr (T::Option == Double) {
            /**
             * Double version is implemented with gamma = 6 and N = 9 [1] to make full use of precision of double
             * 
             * TODO: The implementation is less precise than the implementation in STL, try to find out the difference
             * 
             * Compare:
             * The implementation use different number of series term for float and double, while STL provided only one implementation,
             * this version should run faster.
             */
            constexpr static int count = 9;
            constexpr static double coeffcients[count]{228.9344030404165, -342.8104127892456, 151.3843107005646, -20.01174920149977, 0.4619036553182262, -0.0001214195995667437, -1.535239091824004E-6, 1.102873029688190E-6, -2.202670452322396E-7};
            T temp = s + T(6.5);
            temp -= (s + T(0.5)) * ln(temp);
            T ser(1.000000000000123);
            T copy(s);
            for (int j = 0; j < count; ++j) {
                copy += T(1);
                ser += T(coeffcients[j]) / copy;
            }
            return -temp + T(0.91893853320467274178) + ln(ser / s);
        }
        else {
            /**
             * float version is implemented with gamma = 3 and N = 4 [1] to make full use of precision of float
             */
            constexpr static int count = 4;
            constexpr static float coeffcients[count]{7.6845130, -3.284681, 0.05832037, 0.0001856071};
            T temp = s + T(3.5);
            temp -= (s + T(0.5)) * ln(temp);
            T ser(0.9999998);
            T copy(s);
            for (int j = 0; j < count; ++j) {
                copy += T(1);
                ser += T(coeffcients[j]) / copy;
            }
            return -temp + T(0.91893853320467274178) + ln(ser / s);
        }
    }

    template<ScalarOption Option>
    inline Real<Option> gamma(const Real<Option>& s) {
        return exp(lnGamma(s));
    }

    template<ScalarOption Option>
    inline Real<Option> beta(const Real<Option>& s1, const Real<Option>& s2) {
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

    template<Scalar T>
    T bigamma(const T& x, const T& step) {
        return Differential<T>::ridders(lnGamma<T>, x, step);
    }

    template<ScalarOption Option>
    Real<Option> erf(const Real<Option>& x) {
        using T = Real<Option>;
        T x2 = square(x);
        return (x.isNegative()) ? -gammaP(T(0.5), x2) : gammaP(T(0.5), x2);
    }

    template<Scalar T>
    T erfc(const T& x) {
        T x2 = square(x);
        return (x.isNegative()) ? (T(1) + gammaP(T(0.5), x2)) : gammaQ(T(0.5), x2);
    }

    template<ScalarOption Option>
    Real<Option> standardNormalDistribution(const Real<Option>& x) {
        using T = Real<Option>;
        return (erf(x / sqrt(T(2))) + T(1)) >> 1U;
    }
}
