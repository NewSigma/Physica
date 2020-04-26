#include "../../Header/Integrate.h"

Numerical rectangular(Numerical func(const Numerical&), const Numerical& x0, const Numerical& x1) {
    Numerical deltaY = getZero();

    Numerical point(x0);
    while(point < x1) {
        Numerical temp = func(point);
        deltaY += temp;
        point += BasicConst::getInstance().getStepSize();
    }
    return deltaY * BasicConst::getInstance().getStepSize();
}

Numerical ladder(Numerical func(const Numerical&), const Numerical& x0, const Numerical& x1) {
    Numerical deltaY = (x0 + x1) / BasicConst::getInstance().get_2();
    Numerical point = x0 + BasicConst::getInstance().getStepSize();
    while(point < x1) {
        deltaY += func(point);
        point += BasicConst::getInstance().getStepSize();
    }
    return deltaY * BasicConst::getInstance().getStepSize();
}

Numerical simpson(Numerical func(const Numerical&), const Numerical& x0, const Numerical& x1) {
    Numerical temp0 = x0 + x1, temp1 = getZero(), temp2 = getZero();
    bool odd = true;
    Numerical point = x0 + BasicConst::getInstance().getStepSize();
    while(point < x1) {
        if(odd)
            temp1 += func(point);
        else
            temp2 += func(point);
        odd = !odd;
        point += BasicConst::getInstance().getStepSize();
    }
    return (temp0 + temp1 + BasicConst::getInstance().get_4() + temp2 * BasicConst::getInstance().get_2()) * BasicConst::getInstance().getStepSize() / BasicConst::getInstance().get_3();
}