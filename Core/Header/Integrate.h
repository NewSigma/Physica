#ifndef PHYSICA_INTEGRATE_H
#define PHYSICA_INTEGRATE_H

#include "Function.h"

namespace Physica::Core {
    class Integrate {
        Function func;
        Scalar from;
        Scalar to;
        Scalar stepSize;
    public:
        //Reference: Scalar Recipes in C++
        enum IntegrateMethod {
            Rectangular,
            Ladder,
            Simpson,
            Simpson_3_8,
            Bode
        };
        Integrate(Function func, Scalar from, Scalar to, Scalar stepSize = BasicConst::getInstance().getStepSize());
        [[nodiscard]] Scalar solve(IntegrateMethod method);
    };
}

#endif
