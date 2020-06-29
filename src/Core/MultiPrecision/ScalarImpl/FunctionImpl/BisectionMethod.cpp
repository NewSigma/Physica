/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    template<>
    Scalar<MultiPrecision, false> bisectionMethod(
            Scalar<MultiPrecision, false> func(const Scalar<MultiPrecision, false>&)
            , const Scalar<MultiPrecision, false>& n
            , const Scalar<MultiPrecision, false>& x1, const Scalar<MultiPrecision, false>& x2
            , const Scalar<MultiPrecision, false>& y1, const Scalar<MultiPrecision, false>& y2) {
        if(n == y1)
            return Scalar<MultiPrecision, false>(x1);
        if(n == y2)
            return Scalar<MultiPrecision, false>(x2);

        if(((n - y1).getLength() ^ (n - y2).getLength()) >= 0) // NOLINT(hicpp-signed-bitwise)
            qFatal("Root is nonexistent.");

        Scalar<MultiPrecision, false> result = (x1 + x2) >> 1;
        Scalar<MultiPrecision, false> y_result(static_cast<SignedScalarUnit>(1));

        Scalar<MultiPrecision, false> error = (y1 - y2).toAbs() >> 1;
        Scalar<MultiPrecision, false> x_left(x1);
        Scalar<MultiPrecision, false> x_right(x2);
        Scalar<MultiPrecision, false> y_left(y1);

        bool delta_left_sign = (n - y_left).getLength() > 0;
        bool delta_right_sign;
        do {
            y_result = func(result);
            delta_right_sign = (n - y_result).getLength() > 0;

            if(delta_left_sign == delta_right_sign) {
                x_left = result;
                y_left = y_result;
                delta_left_sign = (n - y_left).getLength() > 0;
            }
            else
                x_right = result;
            result = (x_left + x_right) >> 1;
            error >>= 1;
        } while(result.getPower() - error.getPower() < GlobalPrecision);
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> bisectionMethod(
            Scalar<MultiPrecision, false> func(const Scalar<MultiPrecision, false>&)
            , const Scalar<MultiPrecision, false>& n
            , const Scalar<MultiPrecision, false>& x1, const Scalar<MultiPrecision, false>& x2
            , const Scalar<MultiPrecision, false>& y1, const Scalar<MultiPrecision, false>& y2) {
        if(n == y1)
            return Scalar<MultiPrecision, true>(x1);
        if(n == y2)
            return Scalar<MultiPrecision, true>(x2);

        if(((n - y1).getLength() ^ (n - y2).getLength()) >= 0) // NOLINT(hicpp-signed-bitwise)
            qFatal("Root is nonexistent.");

        Scalar<MultiPrecision, false> result = (x1 + x2) >> 1;
        Scalar<MultiPrecision, false> y_result(static_cast<SignedScalarUnit>(1));

        Scalar<MultiPrecision, false> error = (y1 - y2).toAbs() >> 1;
        Scalar<MultiPrecision, false> x_left(x1);
        Scalar<MultiPrecision, false> x_right(x2);
        Scalar<MultiPrecision, false> y_left(y1);

        bool delta_left_sign = (n - y_left).getLength() > 0;
        bool delta_right_sign;
        do {
            y_result = func(result);
            delta_right_sign = (n - y_result).getLength() > 0;

            if(delta_left_sign == delta_right_sign) {
                x_left = result;
                y_left = y_result;
                delta_left_sign = (n - y_left).getLength() > 0;
            }
            else
                x_right = result;
            result = (x_left + x_right) >> 1;
            error >>= 1;
        } while(result.getPower() - error.getPower() < GlobalPrecision);
        return Scalar<MultiPrecision, true>(result).toUnitA();
    }
}