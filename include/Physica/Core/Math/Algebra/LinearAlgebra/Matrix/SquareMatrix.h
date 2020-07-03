/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SQUAREMATRIX_H
#define PHYSICA_SQUAREMATRIX_H

#include "Matrix.h"

namespace Physica::Core {
    template<class T, MatrixType type, size_t maxSize>
    class SquareMatrix : public Matrix<T, type, maxSize, maxSize> {
    public:
        enum DeterminateMethod {
            GaussMethod,
            LUMethod
        };
    public:
        SquareMatrix();
        explicit SquareMatrix(size_t length);
        SquareMatrix(const SquareMatrix& m) = default;
        SquareMatrix(SquareMatrix&& m) noexcept;
        ~SquareMatrix() = default;
        /* Operators */
        SquareMatrix& operator=(const SquareMatrix& m) = default;
        SquareMatrix& operator=(SquareMatrix&& m) noexcept;
        /* Getters */
        [[nodiscard]] T determinate(DeterminateMethod method);

        static SquareMatrix getUnitMatrix(size_t length);
    };
}

#include "SquareMatrixImpl.h"

#endif