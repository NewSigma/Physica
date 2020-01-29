#ifndef _Physica_C_ClimbMountainAlgorithm_H
#define _Physica_C_ClimbMountainAlgorithm_H

#include "RealNumber.h"

class HillClimbingAlgorithm {
public:
    HillClimbingAlgorithm(RealNumberA* func(RealNumber*), RealNumber* x_initial, RealNumber* stepSize);
	~HillClimbingAlgorithm();
    void getExtremal();
    RealNumber* getMinStep();
	void setMinStep(RealNumber* minStep);
private:
    RealNumberA* (*func)(RealNumber*);
    RealNumber* x_initial;
    RealNumber* stepSize;
    RealNumber* minStep;
};

#endif