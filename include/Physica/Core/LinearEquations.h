/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_LINEAREQUATIONS_H
#define PHYSICA_LINEAREQUATIONS_H

#include "Matrix.h"

namespace Physica::Core {
    class LinearEquations {
    public:
        enum LinearEquationsMethod {
            GaussJordanPartial,
            GaussJordanComplete,
            GaussEliminationPartial,
            GaussEliminationComplete,
            LUMethod
        };
    private:
        Matrix& matrix;
    public:
        explicit LinearEquations(Matrix& m) noexcept;

        const Vector& solve(LinearEquationsMethod method);
    };
}

#endif