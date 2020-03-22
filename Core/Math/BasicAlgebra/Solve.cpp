#include "../../Header/Solve.h"
#include "../../Header/Const.h"
/*
 * Find the numerical root of algebra equation.
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
/*
 * y1 = func(x1), y2 = func(x2)
 * Solve the x when n = func(x)
 * May return nullptr.
 */
Numerical* bisectionMethod(Numerical* func(const Numerical&), const Numerical& n, const Numerical& x1, const Numerical& x2)    {
    auto y1 = func(x1);
    auto y2 = func(x2);
    auto result = bisectionMethod(func, n, x1, x2, *y1, *y2);
    delete y1;
    delete y2;
    return result;
}
//May return nullptr.
Numerical* bisectionMethod(Numerical* func(const Numerical&), const Numerical& n, const Numerical& x1, const Numerical& x2, const Numerical& y1, const Numerical& y2) {
    if(n == y1)
        return new Numerical(x1);
    if(n == y2)
        return new Numerical(x2);

    auto delta_left = subtract(n, y1);
    auto delta_right = subtract(n, y2);
    if(delta_left->sign == delta_right->sign)
        return nullptr;
    delete delta_left;
    delete delta_right;

    auto result = add(x1, x2);
    *result /= *const_1->_2;
    Numerical* y_result;

    auto error = subtract(y1, y2);
    error->sign = true;
    *error << *divide(*error, *const_1->_2);

    auto x_left = new Numerical(x1);
    auto x_right = new Numerical(x2);
    auto y_left = new Numerical(y1);

    delta_left = subtract(n, *y_left);
    bool delta_left_sign = delta_left->sign;
    bool delta_right_sign;
    delete delta_left;

    do {
        y_result = func(*result);
        delta_right = subtract(n, *y_result);
        delta_right_sign = delta_right->sign;
        delete delta_right;

        if(delta_left_sign == delta_right_sign) {
            delete x_left;
            delete y_left;
            x_left = result;
            y_left = y_result;

            delta_left = subtract(n, *y_left);
            delta_left_sign = delta_left->sign;
            delete delta_left;
        }
        else {
            delete x_right;
            delete y_result;
            x_right = result;
        }
        result = add(*x_left, *x_right);
        *result << *divide(*result, *const_1->_2);
        *error << *divide(*error, *const_1->_2);
    } while(result->power - error->power < const_1->GlobalPrecision);
    result->a = 1;

    delete x_left;
    delete x_right;
    delete y_left;
    delete error;

    return result;
}