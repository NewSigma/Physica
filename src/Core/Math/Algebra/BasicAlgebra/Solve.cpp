#include <QtCore/qlogging.h>
#include "Physica/Core/Solve.h"
#include "Physica/Core/Scalar.h"
/*
 * Find the numerical root of algebra equation.
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
namespace Physica::Core {
    /*
     * Solve the x when n = func(x)
     */
    Scalar Solve::bisectionMethod(Scalar func(const Scalar&), const Scalar& n
            , const Scalar& x1, const Scalar& x2) {
        Scalar y1 = func(x1);
        Scalar y2 = func(x2);
        return bisectionMethod(func, n, x1, x2, y1, y2);
    }

    Scalar Solve::bisectionMethod(Scalar func(const Scalar&), const Scalar& n
            , const Scalar& x1, const Scalar& x2, const Scalar& y1, const Scalar& y2) {
        if(n == y1)
            return Scalar(x1);
        if(n == y2)
            return Scalar(x2);

        if((Scalar::sub(n, y1).getLength() ^ Scalar::sub(n, y2).getLength()) >= 0) // NOLINT(hicpp-signed-bitwise)
            qFatal("Root is nonexistent.");

        Scalar result = Scalar::div(Scalar::add(x1, x2), BasicConst::getInstance().get_2());
        Scalar y_result = getZero();

        Scalar error = Scalar::div(Scalar::sub(y1, y2).toAbs(), BasicConst::getInstance().get_2());

        Scalar x_left(x1);
        Scalar x_right(x2);
        Scalar y_left(y1);

        bool delta_left_sign = Scalar::sub(n, y_left).getLength() > 0;
        bool delta_right_sign;
        do {
            y_result = func(result);
            delta_right_sign = Scalar::sub(n, y_result).getLength() > 0;

            if(delta_left_sign == delta_right_sign) {
                x_left = result;
                y_left = y_result;
                delta_left_sign = Scalar::sub(n, y_left).getLength() > 0;
            }
            else
                x_right = result;
            result = Scalar::div(Scalar::add(x_left, x_right), BasicConst::getInstance().get_2());
            error = Scalar::div(error, BasicConst::getInstance().get_2());
        } while(result.getPower() - error.getPower() < GlobalPrecision);
        result.toUnitA();

        return result;
    }
}
