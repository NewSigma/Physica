#ifndef _Physica_C_ClimbMountainAlgorithm_H
#define _Physica_C_ClimbMountainAlgorithm_H

namespace Physica::Core {
    class Numerical;

    class HillClimbingAlgorithm {
    public:
        HillClimbingAlgorithm(Numerical* func(Numerical*), Numerical* x_initial, Numerical* stepSize);
        ~HillClimbingAlgorithm();
        void getExtremal();
        Numerical* getMinStep();
        void setMinStep(Numerical* minStep);
    private:
        Numerical* (*func)(Numerical*);
        Numerical* x_initial;
        Numerical* stepSize;
        Numerical* minStep;
    };
}

#endif