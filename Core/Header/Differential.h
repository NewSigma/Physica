#ifndef PHYSICA_DIFFERENTIAL_H
#define PHYSICA_DIFFERENTIAL_H

#include "FunctionTree.h"

namespace Physica::Core {
    class Numerical;

    class Differential {
        FunctionTree func;
        Numerical at;
        Numerical stepSize;
    public:
        //Reference: Numerical Recipes in C++
        enum DifferentialMethod {
            DoublePoint,
            Forward,
            Backward
        };
        Differential(FunctionTree func, Numerical at, Numerical stepSize = BasicConst::getInstance().getStepSize());
        [[nodiscard]] Numerical solve(DifferentialMethod method) const;
    };
}

#endif
