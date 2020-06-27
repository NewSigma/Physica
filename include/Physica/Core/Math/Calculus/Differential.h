#ifndef PHYSICA_DIFFERENTIAL_H
#define PHYSICA_DIFFERENTIAL_H

#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/Math/Calculus/Function.h"

namespace Physica::Core {
    class Differential {
        Function func;
        MultiScalar at;
        MultiScalar stepSize;
    public:
        //Reference: Scalar Recipes in C++
        enum DifferentialMethod {
            DoublePoint,
            Forward,
            Backward
        };
        Differential(Function func, MultiScalar at
                , MultiScalar stepSize = MultiScalar(BasicConst::getInstance().getStepSize()));
        [[nodiscard]] MultiScalar solve(DifferentialMethod method);
    };
}

#endif
