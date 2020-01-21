#ifndef _Physica_C_Const_H
#define _Physica_C_Const_H

#include "RealNumber.h"

class Const_1 {
public:
    int MachinePrecision;
    const RealNumberA* R_MAX;
    const RealNumberA* ZERO;
    const RealNumberA* ONE;
    const RealNumberA* MINUS_ONE;
    const RealNumberA* TWO;
    const RealNumberA* MINUS_TWO;
    const RealNumberA* PI;

    Const_1();
    ~Const_1();

    RealNumberA* getZero() const { return new RealNumberA(ZERO); }
    RealNumberA* getOne() const { return new RealNumberA(ONE); }
    RealNumberA* getTwo() const { return new RealNumberA(TWO); }
};

class Const_2 {
public:
    const RealNumberA* PI_DIVIDE_TWO;
    const RealNumberA* MINUS_PI_DIVIDE_TWO;

    Const_2();
    ~Const_2();
};

#endif