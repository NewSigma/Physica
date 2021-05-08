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
#ifndef PHYSICA_SPATIALFUNCTIONS_H
#define PHYSICA_SPATIALFUNCTIONS_H

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

namespace Physica::Core {
    /**
     * Return the logarithm of gamma(s). s must be positive.
     * 
     * Implemented with gamma = 6 and N = 9 [1] to make full use of precision of double
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
        return -temp + ln(T(2.5066282746310005) * ser / s);
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
        return -temp + ln(T(2.506628) * ser / s);
    }
    /**
     * Return the logarithm of gamma(s). s must be positive.
     *
     * Reference:
     * [1] Helmut Werner and Robert Collinge.Chebyshev Approximations to the Gamma Function[J]
     *     .Mathematics of Computation, Vol. 15, No. 74 (Apr., 1961), pp. 195-197
     */
    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> lnGamma(const Scalar<MultiPrecision, errorTrack>& s) {
        Q_ASSERT(s.isPositive());
        Scalar<MultiPrecision, errorTrack> result = Scalar<MultiPrecision, errorTrack>::getZero();
        Scalar<MultiPrecision, false> count = Scalar<MultiPrecision, false>::getOne();
        Scalar<MultiPrecision, errorTrack> temp;
        //Handle count = 1
        temp = count << 1;
        Scalar<MultiPrecision, false> doubleCount_1 = temp - Scalar<MultiPrecision, false>(1);
        temp *= doubleCount_1 * (s ^ doubleCount_1);
        temp = bernoulli(count) / temp;
        bool flag = true;
        ++count;
        bool minus = true;
        while(flag) {
            result += minus ? -temp : temp; //Optimize: new helper function Scalar::toOpposite(bool) to avoid branch

            temp = count << 1;
            doubleCount_1 = temp - Scalar<MultiPrecision, false>(1);
            temp *= doubleCount_1 * (s ^ doubleCount_1);
            temp = bernoulli(count) / temp;
            flag = (temp.getPower() - result.getPower()) > GlobalPrecision;
            ++count;
            minus = !minus;
        }
        result.toUnitA();

        const auto term1 = (s - (Scalar<MultiPrecision, errorTrack>::getOne() >> 1)) * ln(s);
        const auto term2 = ln(MathConst::getInstance().PI << 1) >> 1;
        result = term1 - s + term2 + result;
        return result;
    }

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> gamma(const Scalar<type, errorTrack>& s) {
        return exp(lnGamma(s));
    }
}

#endif