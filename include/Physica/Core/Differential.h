#ifndef PHYSICA_DIFFERENTIAL_H
#define PHYSICA_DIFFERENTIAL_H

#include "Scalar.h"
#include "Function.h"

namespace Physica::Core {
    class Differential {
        Function func;
        Scalar at;
        Scalar stepSize;
    public:
        //Reference: Scalar Recipes in C++
        enum DifferentialMethod {
            DoublePoint,
            Forward,
            Backward
        };
        Differential(Function func, Scalar at, Scalar stepSize = BasicConst::getInstance().getStepSize());
        [[nodiscard]] Scalar solve(DifferentialMethod method);
    };
}

#endif
