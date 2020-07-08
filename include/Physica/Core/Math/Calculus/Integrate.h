#ifndef PHYSICA_INTEGRATE_H
#define PHYSICA_INTEGRATE_H

#include "Physica/Core/Math/Calculus/Function/Function.h"
#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    class Integrate {
        Function func;
        MultiScalar from;
        MultiScalar to;
        Scalar<MultiPrecision, false> stepSize;
    public:
        //Reference: MultiScalar Recipes in C++
        enum IntegrateMethod {
            Rectangular,
            Ladder,
            Simpson,
            Simpson_3_8,
            Bode
        };
        Integrate(Function func, MultiScalar from, MultiScalar to
                , Scalar<MultiPrecision, false> stepSize = BasicConst::getInstance().getStepSize());
        [[nodiscard]] MultiScalar solve(IntegrateMethod method);
    };
}

#endif
