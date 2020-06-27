/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_POLYNOMIAL_H
#define PHYSICA_POLYNOMIAL_H

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    class Polynomial {
    private:
        int length;
        const MultiScalar* variable;
        const MultiScalar** coefficients;
    public:
        Polynomial(int length, const MultiScalar* variable);
        ~Polynomial();
        void addCoefficients(int index, const MultiScalar* n);
        MultiScalar calculate();
    };
}

#endif
