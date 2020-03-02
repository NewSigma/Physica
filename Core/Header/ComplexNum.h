/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_COMPLEXNUM_H
#define PHYSICA_COMPLEXNUM_H

#include "Number.h"
#include "RealNumber.h"

class ComplexNum : public Number{
private:
    class RealNumber* real;
    class RealNumber* imagine;
public:
    ComplexNum(class RealNumber* real, class RealNumber* imagine);
    ~ComplexNum();
};

#endif
