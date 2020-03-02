/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_REALINF_H
#define PHYSICA_REALINF_H

#include "Number.h"

class RealInf : public Number {
private:
    bool sign;
public:
    RealInf(bool b);
};

#endif
