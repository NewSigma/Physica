#include "Core/Header/Numerical.h"

namespace Physica::Core {
    Numerical D_double_point(Numerical func(const Numerical&), const Numerical& x0) {
        Numerical x1 = x0 + BasicConst::getInstance().getStepSize();
        Numerical x2 = x0 - BasicConst::getInstance().getStepSize();
        Numerical temp = func(x2) - func(x1);
        return temp / (x2 - x1);
    }

    Numerical D_right(Numerical func(const Numerical&), const Numerical& x0) {
        return (func(x0 + BasicConst::getInstance().getStepSize()) - func(x0)) / BasicConst::getInstance().getStepSize();
    }

    Numerical D_left(Numerical func(const Numerical&), const Numerical& x0) {
        return (func(x0) - func(x0 - BasicConst::getInstance().getStepSize())) / BasicConst::getInstance().getStepSize();
    }
}