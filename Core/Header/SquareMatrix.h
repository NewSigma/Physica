/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SQUAREMATRIX_H
#define PHYSICA_SQUAREMATRIX_H

#include "Matrix.h"

namespace Physica::Core {
    /////////////////////////////////////SquareMatrix//////////////////////////////////////////
    class SquareMatrix : virtual public Matrix {
    public:
        enum SquareMatrixMethod {
            GaussMethod,
            LUMethod
        };
    public:
        void toUnit();
        /* Getters */
        [[nodiscard]] Numerical determinate(SquareMatrixMethod method);
    protected:
        /* Friends */
        friend class LinearEquations;
    };
    /////////////////////////////////////ColumnSquareMatrix//////////////////////////////////////////
    class ColumnSquareMatrix : public ColumnMatrix, public SquareMatrix {
    public:
        ColumnSquareMatrix();
        explicit ColumnSquareMatrix(size_t size);
        ColumnSquareMatrix(const ColumnSquareMatrix& matrix);
        /* Helpers */
        [[nodiscard]] size_t row() const override { return length; }
    };
    /////////////////////////////////////RowSquareMatrix//////////////////////////////////////////
    class RowSquareMatrix : public RowMatrix, public SquareMatrix {
    public:
        RowSquareMatrix();
        explicit RowSquareMatrix(size_t size);
        RowSquareMatrix(const RowSquareMatrix& matrix);
        /* Helpers */
        [[nodiscard]] size_t column() const override { return length; }
    };
}

#endif
