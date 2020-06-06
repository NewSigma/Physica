#ifndef PHYSICA_INTEGRATE_H
#define PHYSICA_INTEGRATE_H

#include "Function.h"

namespace Physica::Core {
    class Integrate {
        Function func;
        Numerical from;
        Numerical to;
        Numerical stepSize;
    public:
        //Reference: Numerical Recipes in C++
        enum IntegrateMethod {
            Rectangular,
            Ladder,
            Simpson,
            Simpson_3_8,
            Bode
        };
        Integrate(Function func, Numerical from, Numerical to, Numerical stepSize = BasicConst::getInstance().getStepSize());
        [[nodiscard]] Numerical solve(IntegrateMethod method);
    };
}

#endif
