#ifndef _Physica_C_ClimbMountainAlgorithm_H
#define _Physica_C_ClimbMountainAlgorithm_H

#include "RealNumber.h"

class ClimbMountainAlgorithm {
public:
    ClimbMountainAlgorithm(RealNumber* func(RealNumber*), RealNumber* x_initial, RealNumber* stepSize);
	~ClimbMountainAlgorithm();
    void getExtremal();
    RealNumber* getMinStep();
	void setMinStep(RealNumber* minStep);
private:
    RealNumber* (*func)(RealNumber*);
    RealNumber* x_initial;
    RealNumber* stepSize;
    RealNumber* minStep;
};

#endif