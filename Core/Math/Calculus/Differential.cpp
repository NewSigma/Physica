#include "../../Header/Differential.h"
#include "../../Header/Const.h"

extern Const_1 const_1;

RealNumber* derivative_right(RealNumber* func(const RealNumber&), const RealNumber& x0) {
    auto x1 = x0 + *const_1.StepSize;
    auto y0 = func(x0);
    auto y1 = func(*x1);
    auto result = *y1 - *y0;
    *result /= *const_1.StepSize;

    delete x1;
    delete y0;
    delete y1;
    return result;
}

RealNumber* derivative_left(RealNumber* func(const RealNumber&), const RealNumber& x0) {
    auto x1 = x0 - *const_1.StepSize;
    auto y0 = func(x0);
    auto y1 = func(*x1);
    auto result = *y0 - *y1;
    *result /= *const_1.StepSize;

    delete x1;
    delete y0;
    delete y1;
    return result;
}