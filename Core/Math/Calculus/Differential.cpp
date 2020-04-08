#include "../../Header/Differential.h"
#include "../../Header/Numerical.h"

Numerical* D_double_point(Numerical* func(const Numerical&), const Numerical& x0) {
    auto x1 = x0 + basicConst->getStepSize();
    auto x2 = x0 - basicConst->getStepSize();
    auto y1 = func(*x1);
    auto y2 = func(*x2);
    auto result = *y2 - *y1;
    auto temp = *x2 - *x1;
    *result /= *temp;
    delete x1;
    delete x2;
    delete y1;
    delete y2;
    delete temp;
    return result;
}

Numerical* D_right(Numerical* func(const Numerical&), const Numerical& x0) {
    auto x1 = x0 + basicConst->getStepSize();
    auto y0 = func(x0);
    auto y1 = func(*x1);
    auto result = *y1 - *y0;
    *result /= basicConst->getStepSize();
    delete x1;
    delete y0;
    delete y1;
    return result;
}

Numerical* D_left(Numerical* func(const Numerical&), const Numerical& x0) {
    auto x1 = x0 - basicConst->getStepSize();
    auto y0 = func(x0);
    auto y1 = func(*x1);
    auto result = *y0 - *y1;
    *result /= basicConst->getStepSize();
    delete x1;
    delete y0;
    delete y1;
    return result;
}