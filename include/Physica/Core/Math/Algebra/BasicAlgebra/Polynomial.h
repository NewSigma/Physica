/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_POLYNOMIAL_H
#define PHYSICA_POLYNOMIAL_H

namespace Physica::Core {
    class Scalar;

    class Polynomial {
    private:
        int length;
        const Scalar* variable;
        const Scalar** coefficients;
    public:
        Polynomial(int length, const Scalar* variable);
        ~Polynomial();
        void addCoefficients(int index, const Scalar* n);
        Scalar* calculate();
    };
}

#endif
