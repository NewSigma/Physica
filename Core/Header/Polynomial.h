/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_POLYNOMIAL_H
#define PHYSICA_POLYNOMIAL_H

#include "Numerical.h"

class Polynomial {
private:
    int length;
    const Numerical* variable;
    const Numerical** coefficients;
public:
    Polynomial(int length, const Numerical* variable);
    ~Polynomial();
    void addCoefficients(int index, const Numerical* n);
    Numerical* calculate();
};

#endif
