#include <QtCore/qlogging.h>
#include "Solve.h"
#include "Numerical.h"
/*
 * Find the numerical root of algebra equation.
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
/*
 * Solve the x when n = func(x)
 */
Numerical bisectionMethod(Numerical func(const Numerical&), const Numerical& n, const Numerical& x1, const Numerical& x2)    {
    Numerical y1 = func(x1);
    Numerical y2 = func(x2);
    return bisectionMethod(func, n, x1, x2, y1, y2);
}

Numerical bisectionMethod(Numerical func(const Numerical&), const Numerical& n, const Numerical& x1, const Numerical& x2, const Numerical& y1, const Numerical& y2) {
    if(n == y1)
        return Numerical(x1);
    if(n == y2)
        return Numerical(x2);

    if((sub(n, y1).getLength() ^ sub(n, y2).getLength()) > 0) // NOLINT(hicpp-signed-bitwise)
        qFatal("Root is nonexistent.");

    Numerical result = add(x1, x2);
    result /= basicConst->get_2();
    Numerical y_result = getZero();

    Numerical error = div(sub(y1, y2).toAbs(), basicConst->get_2());

    Numerical x_left(x1);
    Numerical x_right(x2);
    Numerical y_left(y1);

    bool delta_left_sign = sub(n, y_left).getLength() > 0;
    bool delta_right_sign;
    do {
        y_result = func(result);
        delta_right_sign = sub(n, y_result).getLength() > 0;

        if(delta_left_sign == delta_right_sign) {
            x_left = result;
            y_left = y_result;
            delta_left_sign = sub(n, y_left).getLength() > 0;
        }
        else
            x_right = result;
        result = div(add(x_left, x_right), basicConst->get_2());
        error = div(error, basicConst->get_2());
    } while(result.getPower() - error.getPower() < basicConst->getGlobalPrecision());
    result.toUnitA();

    return result;
}