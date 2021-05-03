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
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.32
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> lnGamma(const Scalar<type, errorTrack>& s) {
        using T = Scalar<type, false>;
        assert(s.isPositive());
        constexpr static int count = 6;
        constexpr static double coeffcients[count]{76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179E-2, -0.5395239384953E-5};
        Scalar<type, errorTrack> temp = s + T(5.5);
        temp -= (s + T(0.5)) * ln(temp);
        Scalar<type, errorTrack> ser(1.000000000190015);
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