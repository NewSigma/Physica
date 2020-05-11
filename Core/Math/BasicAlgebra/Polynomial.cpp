/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <iostream>
#include "Core/Header/Polynomial.h"
#include "Core/Header/Numerical.h"

namespace Physica::Core {
    /*
     * A polynomial equals to a0 + a1 * x + a2 * x^2 + ... + an * x^n.
     * Here length = n, variable = x, coefficients is a array of a1, a2, ..., an.
     */
    Polynomial::Polynomial(int len, const Numerical* var) {
        length = len;
        variable = var;
        coefficients = new const Numerical*[length];
    }

    Polynomial::~Polynomial() {
        delete[] coefficients;
    }

    void Polynomial::addCoefficients(int index, const Numerical* n) {
        if(index > -1 && index < length)
            coefficients[index] = n;
        else
            std::cout << "[Polynomial] Invalid coefficient index." << std::endl;
    }

    Numerical* Polynomial::calculate() {
        auto result = new Numerical(getZero());
        for(int i = length - 1; i > 0; --i) {
            *result += *coefficients[i];
            *result *= *variable;
        }
        *result += *coefficients[0];
        return result;
    }
}
