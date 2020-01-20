/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_POLYNOMIAL_H
#define PHYSICA_POLYNOMIAL_H

#include "RealNumber.h"

class Polynomial {
private:
    int length;
    const RealNumber* variable;
    const RealNumber** coefficients;
public:
    Polynomial(int length, const RealNumber* variable);
    ~Polynomial();
    void addCoefficients(int index, const RealNumber* n);
    RealNumber* calculate();
};

#endif
