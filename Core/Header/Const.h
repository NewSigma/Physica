#ifndef _Physica_C_Const_H
#define _Physica_C_Const_H

#include "RealNum.h"

class Numerical;

class Const_1 {
public:
    //_1 stands by integer 1.
    int GlobalPrecision;
    const Numerical* StepSize;
    const Numerical* R_MAX;
    const Numerical* _0;
    const Numerical* _1;
    const Numerical* Minus_1;
    const Numerical* _2;
    const Numerical* Minus_2;
    const Numerical* _3;
    const Numerical* Minus_3;
    const Numerical* _4;
    const Numerical* Minus_4;

    Const_1();
    ~Const_1();
};

class Const_2 {
public:
    const RealNum* _0;
    const RealNum* _1;
    const Numerical* PI;
    const Numerical* E;
    //Here PI_2 stands by PI / 2.
    const Numerical* PI_2;
    const Numerical* Minus_PI_2;

    Const_2();
    ~Const_2();
};

#endif