#ifndef _Physica_C_Const_H
#define _Physica_C_Const_H

#include "RealNumber.h"

class Const_1 {
public:
    //_1 stands by integer 1.
    int GlobalPrecision;
    const RealNumber* StepSize;
    const RealNumber* R_MAX;
    const RealNumber* _0;
    const RealNumber* _1;
    const RealNumber* Minus_1;
    const RealNumber* _2;
    const RealNumber* Minus_2;
    const RealNumber* _3;
    const RealNumber* Minus_3;
    const RealNumber* _4;
    const RealNumber* Minus_4;

    Const_1();
    ~Const_1();

    RealNumber* getZero() const { return new RealNumber(_0); }
    RealNumber* getOne() const { return new RealNumber(_1); }
    RealNumber* getTwo() const { return new RealNumber(_2); }
};

class Const_2 {
public:
    const RealNumber* PI;
    const RealNumber* E;
    //Here PI_2 stands by PI / 2.
    const RealNumber* PI_2;
    const RealNumber* Minus_PI_2;

    Const_2();
    ~Const_2();
};

#endif