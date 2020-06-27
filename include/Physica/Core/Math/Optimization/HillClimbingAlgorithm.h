#ifndef _Physica_C_ClimbMountainAlgorithm_H
#define _Physica_C_ClimbMountainAlgorithm_H

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    class HillClimbingAlgorithm {
    public:
        HillClimbingAlgorithm(MultiScalar* func(MultiScalar*), MultiScalar* x_initial, MultiScalar* stepSize);
        ~HillClimbingAlgorithm();
        void getExtremal();
        MultiScalar* getMinStep();
        void setMinStep(MultiScalar* minStep);
    private:
        MultiScalar* (*func)(MultiScalar*);
        MultiScalar* x_initial;
        MultiScalar* stepSize;
        MultiScalar* minStep;
    };
}

#endif