#include "../../Header/Differential.h"
#include "../../Header/Numerical.h"

Numerical D_double_point(Numerical func(const Numerical&), const Numerical& x0) {
    Numerical x1 = x0 + basicConst->getStepSize();
    Numerical x2 = x0 - basicConst->getStepSize();
    Numerical temp = func(x2) - func(x1);
    return temp / (x2 - x1);
}

Numerical D_right(Numerical func(const Numerical&), const Numerical& x0) {
    return (func(x0 + basicConst->getStepSize()) - func(x0)) / basicConst->getStepSize();
}

Numerical D_left(Numerical func(const Numerical&), const Numerical& x0) {
    return (func(x0) - func(x0 - basicConst->getStepSize())) / basicConst->getStepSize();
}