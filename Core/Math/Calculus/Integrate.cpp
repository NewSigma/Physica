#include "../../Header/Integrate.h"

Numerical rectangular(Numerical func(const Numerical&), const Numerical& x0, const Numerical& x1) {
    Numerical deltaY = getZero();

    Numerical point(x0);
    while(point < x1) {
        Numerical temp = func(point);
        deltaY += temp;
        point += basicConst->getStepSize();
    }
    return deltaY * basicConst->getStepSize();
}

Numerical ladder(Numerical func(const Numerical&), const Numerical& x0, const Numerical& x1) {
    Numerical deltaY = (x0 + x1) / basicConst->get_2();
    Numerical point = x0 + basicConst->getStepSize();
    while(point < x1) {
        deltaY += func(point);
        point += basicConst->getStepSize();
    }
    return deltaY * basicConst->getStepSize();
}

Numerical simpson(Numerical func(const Numerical&), const Numerical& x0, const Numerical& x1) {
    Numerical temp0 = x0 + x1, temp1 = getZero(), temp2 = getZero();
    bool odd = true;
    Numerical point = x0 + basicConst->getStepSize();
    while(point < x1) {
        if(odd)
            temp1 += func(point);
        else
            temp2 += func(point);
        odd = !odd;
        point += basicConst->getStepSize();
    }
    return (temp0 + temp1 + basicConst->get_4() + temp2 * basicConst->get_2()) * basicConst->getStepSize() / basicConst->get_3();
}