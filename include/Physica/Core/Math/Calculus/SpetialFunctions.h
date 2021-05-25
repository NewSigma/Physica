/*
 * Copyright 2020-2021 WeiBo He.
 *
 * This file is part of Physica.

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
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

namespace Physica::Core {
    /**
     * Return the logarithm of gamma(s). s must be positive.
     * 
     * Implemented with gamma = 6 and N = 9 [1] to make full use of precision of double
     * 
     * TODO: The implementation is less precise than the implementation in STL, try to find out the difference
     * 
     * Compare:
     * The implementation use different number of series term for float and double, while STL provided only one implementation,
     * this version should run faster.
     * 
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.156
     * [2] Lanczos, C. 1964, SIAM Journal on Numerical Analysis, ser. B, vol. 1, pp. 86-96
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> lnGamma(const Scalar<type, errorTrack>& s) {
        using T = Scalar<type, false>;
        assert(s.isPositive());
        constexpr static int count = 9;
        constexpr static double coeffcients[count]{228.9344030404165, -342.8104127892456, 151.3843107005646, -20.01174920149977, 0.4619036553182262, -0.0001214195995667437, -1.535239091824004E-6, 1.102873029688190E-6, -2.202670452322396E-7};
        Scalar<type, errorTrack> temp = s + T(6.5);
        temp -= (s + T(0.5)) * ln(temp);
        Scalar<type, errorTrack> ser(1.000000000000123);
        Scalar<type, errorTrack> copy(s);
        for (int j = 0; j < count; ++j) {
            copy += T(1);
            ser += T(coeffcients[j]) / copy;
        }
        return -temp + T(0.91893853320467274178) + ln(ser / s);
    }
    /**
     * Return the logarithm of gamma(s). s must be positive.
     * 
     * Implemented with gamma = 3 and N = 4 [1] to make full use of precision of float
     *
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.32
     * [2] Lanczos, C. 1964, SIAM Journal on Numerical Analysis, ser. B, vol. 1, pp. 86-96
     */
    template<bool errorTrack>
    Scalar<Float, errorTrack> lnGamma(const Scalar<Float, errorTrack>& s) {
        using T = Scalar<Float, false>;
        assert(s.isPositive());
        constexpr static int count = 4;
        constexpr static float coeffcients[count]{7.6845130, -3.284681, 0.05832037, 0.0001856071};
        Scalar<Float, errorTrack> temp = s + T(3.5);
        temp -= (s + T(0.5)) * ln(temp);
        Scalar<Float, errorTrack> ser(0.9999998);
        Scalar<Float, errorTrack> copy(s);
        for (int j = 0; j < count; ++j) {
            copy += T(1);
            ser += T(coeffcients[j]) / copy;
        }
        return -temp + T(0.91893853320467274178) + ln(ser / s);
    }

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> gamma(const Scalar<type, errorTrack>& s) {
        return exp(lnGamma(s));
    }

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> beta(const Scalar<type, errorTrack>& s1, const Scalar<type, errorTrack>& s2) {
        return exp(lnGamma(s1) + lnGamma(s2) - lnGamma(s1 + s2));
    }

    namespace Internal {
        /**
         * Reference:
         * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
         * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.156
         */
        template<ScalarType type, bool errorTrack>
        Scalar<type, errorTrack> incompGamma1(const Scalar<type, errorTrack>& a, const Scalar<type, errorTrack>& x) {
            using T = Scalar<type, errorTrack>;
            assert(a.isPositive() && !x.isNegative());
            assert(x < a + T(1)); //When x > a + 1, the algorithm is slow, use the other method is better
            T ap = a;
            T temp = reciprocal(a);
            T sum = temp;
            do {
                ap += T(1);
                temp *= x / ap;
                sum += temp;
            } while (fabs(temp.getTrivial()) >= fabs(sum.getTrivial()) * std::numeric_limits<typename T::TrivialType>::epsilon());
            return sum * exp(T(-x + a * ln(x) - lnGamma(a)));
        }
        /**
         * Reference:
         * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
         * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.161
         */
        template<ScalarType type, bool errorTrack>
        Scalar<type, errorTrack> incompGamma2(const Scalar<type, errorTrack>& a, const Scalar<type, errorTrack>& x) {
            using T = Scalar<type, errorTrack>;
            assert(a.isPositive() && !x.isNegative());
            assert(x > a + T(1)); //When x < a + 1, the algorithm is slow, use the other method is better
            constexpr typename T::TrivialType epsilon = std::numeric_limits<typename T::TrivialType>::epsilon();
            constexpr typename T::TrivialType floatMin = std::numeric_limits<typename T::TrivialType>::min() / epsilon;

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
                if (copy_d.toAbs().getTrivial() < floatMin)
                    d = T(floatMin);

                c = an / c + b;
                T copy_c(c);
                if (copy_c.toAbs().getTrivial() < floatMin)
                    c = T(floatMin);
                d = reciprocal(d);
                temp = c * d;
                h *= temp; 
            } while (T(temp - T(1)).toAbs().getTrivial() >= epsilon);
            return h * exp(T(-x + a * ln(x) - lnGamma(a)));
        }
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> gammaP(const Scalar<type, errorTrack>& a, const Scalar<type, errorTrack>& x) {
        assert(a.isPositive() && !x.isNegative());
        using T = Scalar<type, errorTrack>;
        return (x < a + T(1)) ? Internal::incompGamma1(a, x) : (T(1) - Internal::incompGamma2(a, x));
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> gammaQ(const Scalar<type, errorTrack>& a, const Scalar<type, errorTrack>& x) {
        assert(a.isPositive() && !x.isNegative());
        using T = Scalar<type, errorTrack>;
        return (x < a + T(1)) ? (T(1) - Internal::incompGamma1(a, x)) : Internal::incompGamma2(a, x);
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> erf(const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T x2 = square(x);
        return (x.isNegative()) ? -gammaP(T(0.5), x2) : gammaP(T(0.5), x2);
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> erfc(const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T x2 = square(x);
        return (x.isNegative()) ? (T(1) + gammaP(T(0.5), x2)) : gammaQ(T(0.5), x2);
    }
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.171
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJ0(const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T ax = x;
        ax.toAbs();
        T y, ans1, ans2;
        const T eight = T(8);
        if (ax < eight) {
            y = square(ax);
            ans1 = T(57568490574.0) + y * (T(-13362590354.0) + y * (T(651619640.7) + y * (T(-11214424.18) + y * (T(77392.33017) + y * T(-184.9052456)))));
            ans2 = T(57568490411.0) + y * (T(1029532985.0) + y * (T(9494680.718) + y * (T(59272.64853) + y * (T(267.8532712) + y))));
            return ans1 / ans2;
        }
        else {
            const T z = eight / ax;
            y = square(z);
            ans1 = T(1) + y * (T(-0.1098628627E-2) + y * (T(0.2734510407E-4) + y * (T(-0.2073370639E-5) + y * T(0.2093887211E-6))));
            ans2 = T(-0.1562499995E-1) + y * (T(0.1430488765E-3) + y * (T(-0.6911147651E-5) + y * (T(0.7621095161E-6) - y * T(0.934945152E-7))));
            const T xx = ax - T(0.785398164);
            return sqrt(T(0.636619772) / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
        }
    }
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.171
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJ1(const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        T ax = x;
        ax.toAbs();
        T y, ans1, ans2;
        const T eight = T(8);
        if (ax < eight) {
            y = square(ax);
            ans1 = x * (T(72362614232.0) + y * (T(-7895059235.0) + y * (T(242396853.1) + y * (T(-2972611.439) + y * (T(15704.48260) + y * T(-30.16036606))))));
            ans2 = T(144725228442.0) + y * (T(2300535178.0) + y * (T(18583304.74) + y * (T(99447.43394) + y * (T(376.9991397) + y))));
            return ans1 / ans2;
        }
        else {
            const T z = eight / ax;
            y = square(z);
            ans1 = T(1) + y * (T(0.183105E-2) + y * (T(-0.3516396496E-4) + y * (T(0.2457520174E-5) + y * T(-0.240337019E-6))));
            ans2 = T(0.04687499995) + y * (T(-0.2002690873E-3) + y * (T(0.8449199096E-5) + y * (T(-0.88228987E-6) + y * T(0.105787412E-6))));
            const T xx = ax - T(2.356194491);
            const T result = sqrt(T(0.636619772) / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
            return x.isNegative() ? -result : result;
        }
    }
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.173
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJn(const Integer& n, const Scalar<type, errorTrack>& x) {
        using T = Scalar<type, errorTrack>;
        constexpr int iexp = std::numeric_limits<typename T::TrivialType>::max_exponent / 2;
        constexpr int acc = 160;

        assert(n > 1);
        const T ax = T(x).toAbs();
        const T square_ax = square(ax);
        if (square_ax < T(8 * std::numeric_limits<typename T::TrivialType>::min()))
            return T(0);
        
        const T two_x = T(2) / ax;
        if (ax > T(n)) {
            T bjm = besselJ0(ax); //Possible to optimize: ax is possitive
            T result = besselJ1(ax); //Possible to optimize: ax is possitive
            for (Integer i = 1; i < n; ++i) {
                const T temp = T(i) * two_x * result - bjm;
                bjm = std::move(result);
                result = std::move(temp);
            }
            return result;
        }
        else {
            bool do_sum = false;
            Integer i = n + Integer(sqrt(T(n * acc)));
            i.setByte(0, i.getByte()[0] & (MPUnitMax - 1));
            T bjp(0);
            T bj(1);
            T result(0);
            T sum(0);
            for (; i.isPositive(); --i) {
                T temp = T(i) * two_x * bj - bjp;
                bjp = std::move(bj);
                bj = std::move(temp);
                int k;
                std::frexp(bj.getTrivial(), &k);
                if (k > iexp) {
                    bj = T(std::ldexp(bj.getTrivial(), -iexp));
                    bjp = T(std::ldexp(bjp.getTrivial(), -iexp));
                    result = T(std::ldexp(result.getTrivial(), -iexp));
                    sum = T(std::ldexp(sum.getTrivial(), -iexp));
                }
                if (do_sum) //Optimize: extract loop will avoid branch
                    sum += bj;
                do_sum = !do_sum;
                if (i == n)
                    result = bjp;
            }
            sum = T(2) * sum - bj;
            result /= sum;
            return x.isNegative() && i.isOdd() ? -result : result;
        }
    }
}
