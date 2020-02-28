#include "../../Header/Differential.h"
#include "../../Header/Const.h"

extern Const_1 const_1;

RealNumber* D_double_point(RealNumber* func(const RealNumber&), const RealNumber& x0) {
    auto x1 = x0 + *const_1.StepSize;
    auto y1 = func(*x1);
    auto x2 = x0 - *const_1.StepSize;
    auto y2 = func(*x2);
    auto result = *y1 - *y2;
    auto derivative = *const_1.TWO * *const_1.StepSize;
    *result /= *derivative;
    delete x1;
    delete y1;
    delete x2;
    delete y2;
    delete derivative;
    return result;
}

RealNumber* D_right(RealNumber* func(const RealNumber&), const RealNumber& x0) {
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

RealNumber* D_left(RealNumber* func(const RealNumber&), const RealNumber& x0) {
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