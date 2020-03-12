#include "../../Header/Differential.h"
#include "../../Header/Const.h"
#include "../../Header/RealNum.h"
#include "../../Header/Numerical.h"

AbstractNum* D_double_point(AbstractNum* func(const AbstractNum&), const AbstractNum& x0) {
    auto step = new RealNum(new Numerical(const_1->StepSize));
    AbstractNum* result1 = nullptr, *temp = nullptr;
    //Handle result1 and temp alternatively.
    bool handle = true;
    do {
        delete result1;
        delete temp;
        auto x1 = x0 + *step;
        auto x2 = x0 - *step;
        auto y1 = func(*x1);
        auto y2 = func(*x2);
        auto numerator = *y2 - *y1;
        auto denominator = *x2 - *x1;
        if(handle)
            result1 = *numerator / *denominator;
        else
            temp = *numerator / *denominator;

        delete x1;
        delete x2;
        delete y1;
        delete y2;
        delete numerator;
        delete denominator;

        --step->real->power;
        handle = !handle;
    } while(!handle || *result1 != *temp);

    if(result1->getType() == AbstractNum::RealNumber)
        ((RealNum*)result1)->real->applyError(((RealNum*)temp)->real);
    delete step;
    delete temp;
    return result1;
}

AbstractNum* D_right(AbstractNum* func(const AbstractNum&), const AbstractNum& x0) {
    auto step = new RealNum(new Numerical(const_1->StepSize));
    auto y0 = func(x0);
    AbstractNum* result1 = nullptr, *temp = nullptr;
    //Handle result1 and temp alternatively.
    bool handle = true;
    do {
        delete result1;
        delete temp;
        auto x1 = x0 + *step;
        auto y1 = func(*x1);
        auto numerator = *y1 - *y0;
        auto denominator = *x1 - x0;
        if(handle)
            result1 = *numerator / *denominator;
        else
            temp = *numerator / *denominator;

        delete x1;
        delete y1;
        delete numerator;
        delete denominator;

        --step->real->power;
        handle = !handle;
    } while(!handle || *result1 != *temp);

    if(result1->getType() == AbstractNum::RealNumber)
        ((RealNum*)result1)->real->applyError(((RealNum*)temp)->real);
    delete step;
    delete y0;
    delete temp;
    return result1;
}

AbstractNum* D_left(AbstractNum* func(const AbstractNum&), const AbstractNum& x0) {
    auto step = new RealNum(new Numerical(const_1->StepSize));
    auto y0 = func(x0);
    AbstractNum* result1 = nullptr, *temp = nullptr;
    //Handle result1 and temp alternatively.
    bool handle = true;
    do {
        delete result1;
        delete temp;
        auto x1 = x0 - *step;
        auto y1 = func(*x1);
        auto numerator = *y1 - *y0;
        auto denominator = *x1 - x0;
        if(handle)
            result1 = *numerator / *denominator;
        else
            temp = *numerator / *denominator;

        delete x1;
        delete y1;
        delete numerator;
        delete denominator;

        --step->real->power;
        handle = !handle;
    } while(!handle || *result1 != *temp);

    if(result1->getType() == AbstractNum::RealNumber)
        ((RealNum*)result1)->real->applyError(((RealNum*)temp)->real);
    delete step;
    delete y0;
    delete temp;
    return result1;
}