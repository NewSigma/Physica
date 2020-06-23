#ifndef _Physica_C_ClimbMountainAlgorithm_H
#define _Physica_C_ClimbMountainAlgorithm_H

namespace Physica::Core {
    class Scalar;

    class HillClimbingAlgorithm {
    public:
        HillClimbingAlgorithm(Scalar* func(Scalar*), Scalar* x_initial, Scalar* stepSize);
        ~HillClimbingAlgorithm();
        void getExtremal();
        Scalar* getMinStep();
        void setMinStep(Scalar* minStep);
    private:
        Scalar* (*func)(Scalar*);
        Scalar* x_initial;
        Scalar* stepSize;
        Scalar* minStep;
    };
}

#endif