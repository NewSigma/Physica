#ifndef PHYSICA_DIFFERENTIAL_H
#define PHYSICA_DIFFERENTIAL_H

#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunction.h"

namespace Physica::Core {
    class Differential {
        TreeFunction<> func;
        MultiScalar at;
        MultiScalar stepSize;
    public:
        //Reference: Scalar Recipes in C++
        enum DifferentialMethod {
            DoublePoint,
            Forward,
            Backward
        };
        Differential(TreeFunction<> func, MultiScalar at
                , MultiScalar stepSize = MultiScalar(BasicConst::getInstance().getStepSize()));
        [[nodiscard]] MultiScalar solve(DifferentialMethod method);
    };
}

#endif
