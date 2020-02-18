#ifndef _Physica_C_Const_H
#define _Physica_C_Const_H

#include "RealNumber.h"

class Const_1 {
public:
    int GlobalPrecision;
    const RealNumber* StepSize;
    const RealNumber* R_MAX;
    const RealNumber* ZERO;
    const RealNumber* ONE;
    const RealNumber* MINUS_ONE;
    const RealNumber* TWO;
    const RealNumber* MINUS_TWO;
    const RealNumber* THREE;
    const RealNumber* MINUS_THREE;
    const RealNumber* FOUR;
    const RealNumber* MINUS_FOUR;

    Const_1();
    ~Const_1();

    RealNumber* getZero() const { return new RealNumber(ZERO); }
    RealNumber* getOne() const { return new RealNumber(ONE); }
    RealNumber* getTwo() const { return new RealNumber(TWO); }
};

class Const_2 {
public:
    const RealNumber* PI;
    const RealNumber* E;
    const RealNumber* PI_DIVIDE_TWO;
    const RealNumber* MINUS_PI_DIVIDE_TWO;

    Const_2();
    ~Const_2();
};

#endif