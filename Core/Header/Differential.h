#ifndef PHYSICA_DIFFERENTIAL_H
#define PHYSICA_DIFFERENTIAL_H

#include "Numerical.h"
#include "Function.h"

namespace Physica::Core {
    class Differential {
        Function func;
        Numerical at;
        Numerical stepSize;
    public:
        //Reference: Numerical Recipes in C++
        enum DifferentialMethod {
            DoublePoint,
            Forward,
            Backward
        };
        Differential(Function func, Numerical at, Numerical stepSize = BasicConst::getInstance().getStepSize());
        [[nodiscard]] Numerical solve(DifferentialMethod method);
    };
}

#endif
