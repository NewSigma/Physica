/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SPATIALFUNCTIONS_H
#define PHYSICA_SPATIALFUNCTIONS_H

#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

namespace Physica::Core {
    /*!
     * Return the logarithm of gamma(s). s must be positive.
     *
     * Reference:
     * [1] Helmut Werner and Robert Collinge.Chebyshev Approximations to the Gamma Function[J]
     *     .Mathematics of Computation, Vol. 15, No. 74 (Apr., 1961), pp. 195-197
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> lnGamma(const Scalar<type, errorTrack>& s) {
        Q_ASSERT(s.isPositive());
        Scalar<type, errorTrack> result = Scalar<type, errorTrack>::getZero();
        Scalar<type, false> count = Scalar<type, false>::getOne();
        Scalar<type, errorTrack> temp;
        //Handle count = 1
        temp = count << 1;
        Scalar<type, false> doubleCount_1 = temp - 1;
        temp *= doubleCount_1 * (s ^ doubleCount_1);
        temp = bernoulli(count) / temp;
        bool flag = true;
        ++count;
        bool minus = true;
        while(flag) {
            result += minus ? -temp : temp; //Optimize: new helper function Scalar::toOpposite(bool) to avoid branch

            temp = count << 1;
            doubleCount_1 = temp - 1;
            temp *= doubleCount_1 * (s ^ doubleCount_1);
            temp = bernoulli(count) / temp;
            flag = (temp.getTrivial() / result.getTrivial()) < RelativeError;
            ++count;
            minus = !minus;
        }
        result.setA(temp.getTrivial());

        const auto term1 = (s - (Scalar<type, errorTrack>::getOne() >> 1)) * ln(s);
        const auto term2 = ln(MathConst::getInstance().getPI() << 1) >> 1;
        result = term1 - s + term2 + result;
        return result;
    }

    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> lnGamma(const Scalar<MultiPrecision, errorTrack>& s) {
        Q_ASSERT(s.isPositive());
        Scalar<MultiPrecision, errorTrack> result = Scalar<MultiPrecision, errorTrack>::getZero();
        Scalar<MultiPrecision, false> count = Scalar<MultiPrecision, false>::getOne();
        Scalar<MultiPrecision, errorTrack> temp;
        //Handle count = 1
        temp = count << 1;
        Scalar<MultiPrecision, false> doubleCount_1 = temp - 1;
        temp *= doubleCount_1 * (s ^ doubleCount_1);
        temp = bernoulli(count) / temp;
        bool flag = true;
        ++count;
        bool minus = true;
        while(flag) {
            result += minus ? -temp : temp; //Optimize: new helper function Scalar::toOpposite(bool) to avoid branch

            temp = count << 1;
            doubleCount_1 = temp - 1;
            temp *= doubleCount_1 * (s ^ doubleCount_1);
            temp = bernoulli(count) / temp;
            flag = temp.getPower() - result.getPower() > GlobalPrecision;
            ++count;
            minus = !minus;
        }
        result.toUnitA();

        const auto term1 = (s - (Scalar<MultiPrecision, errorTrack>::getOne() >> 1)) * ln(s);
        const auto term2 = ln(MathConst::getInstance().getPI() << 1) >> 1;
        result = term1 - s + term2 + result;
        return result;
    }

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> gamma(const Scalar<type, errorTrack>& s) {
        return exp(lnGamma(s));
    }
}

#endif