/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_COMPLEXNUM_H
#define PHYSICA_COMPLEXNUM_H

#include "Number.h"
#include "RealNumber.h"

class ComplexNum : public Number{
private:
    RealNumber* real;
    RealNumber* imagine;
public:
    ComplexNum(RealNumber* real, RealNumber* imagine);
    ~ComplexNum();
};

#endif
