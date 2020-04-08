#include "../../Header/Integrate.h"

Numerical* rectangular(Numerical* func(const Numerical&), const Numerical& x0, const Numerical& x1) {
    auto result = getZero();

    auto point = new Numerical(x0);
    while(*point < x1) {
        auto temp = func(*point);
        *result += *temp;
        *point += basicConst->getStepSize();
        delete temp;
    }
    delete point;

    *result *= basicConst->getStepSize();
    return result;
}

Numerical* ladder(Numerical* func(const Numerical&), const Numerical& x0, const Numerical& x1) {
    auto result = getZero();
    *result += x0;
    *result += x1;
    *result /= basicConst->get_2();

    auto point = x0 + basicConst->getStepSize();
    while(*point < x1) {
        auto temp = func(*point);
        *result += *temp;
        *point += basicConst->getStepSize();
        delete temp;
    }
    delete point;

    *result *= basicConst->getStepSize();
    return result;
}

Numerical* simpson(Numerical* func(const Numerical&), const Numerical& x0, const Numerical& x1) {
    auto result = getZero(), temp1 = getZero(), temp2 = getZero();
    *result += x0;
    *result += x1;

    bool odd = true;
    auto point = x0 + basicConst->getStepSize();
    while(*point < x1) {
        auto temp = func(*point);
        if(odd)
            *temp1 += *temp;
        else
            *temp2 += *temp;
        odd = !odd;
        *point += basicConst->getStepSize();
        delete temp;
    }
    delete point;

    *temp1 += basicConst->get_4();
    *temp2 *= basicConst->get_2();
    *result += *temp1;
    *result += *temp2;
    *result *= basicConst->getStepSize();
    *result /= basicConst->get_3();
    delete temp1;
    delete temp2;
    return result;
}