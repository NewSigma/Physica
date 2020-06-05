#ifndef PHYSICA_INTEGRATE_H
#define PHYSICA_INTEGRATE_H

#include "FunctionTree.h"

namespace Physica::Core {
    class Numerical;

    class Integrate {
        FunctionTree func;
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
        Integrate(FunctionTree func, Numerical from, Numerical to, Numerical stepSize = BasicConst::getInstance().getStepSize());
        [[nodiscard]] Numerical solve(IntegrateMethod method) const;
    };
}

#endif
