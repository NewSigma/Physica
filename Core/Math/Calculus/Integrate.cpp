#include "../../Header/Integrate.h"

Numerical* rectangular(Numerical* func(Numerical*), Numerical* x0, Numerical* x1) {
    auto result = getZero();

    auto point = new Numerical(x0);
    while(*point < *x1) {
        auto temp = func(point);
        *result += *temp;
        *point += *const_1->stepSize;
        delete temp;
    }
    delete point;

    *result *= *const_1->stepSize;
    return result;
}

Numerical* ladder(Numerical* func(Numerical*), Numerical* x0, Numerical* x1) {
    auto result = getZero();
    *result += *x0;
    *result += *x1;
    *result /= *const_1->_2;

    auto point = *x0 + *const_1->stepSize;
    while(*point < *x1) {
        auto temp = func(point);
        *result += *temp;
        *point += *const_1->stepSize;
        delete temp;
    }
    delete point;

    *result *= *const_1->stepSize;
    return result;
}

Numerical* Simpson(Numerical* func(Numerical*), Numerical* x0, Numerical* x1) {
    auto result = getZero(), temp1 = getZero(), temp2 = getZero();
    *result += *x0;
    *result += *x1;

    bool odd = true;
    auto point = *x0 + *const_1->stepSize;
    while(*point < *x1) {
        auto temp = func(point);
        if(odd)
            *temp1 += *temp;
        else
            *temp2 += *temp;
        odd = !odd;
        *point += *const_1->stepSize;
        delete temp;
    }
    delete point;

    *temp1 += *const_1->_4;
    *temp2 *= *const_1->_2;
    *result += *temp1;
    *result += *temp2;
    *result *= *const_1->stepSize;
    *result /= *const_1->_3;
    delete temp1;
    delete temp2;
    return result;
}