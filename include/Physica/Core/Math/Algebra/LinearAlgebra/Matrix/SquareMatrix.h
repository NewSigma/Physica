/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SQUAREMATRIX_H
#define PHYSICA_SQUAREMATRIX_H

#include "ColumnMatrix.h"
#include "RowMatrix.h"

namespace Physica::Core {
    /////////////////////////////////////SquareMatrix//////////////////////////////////////////
    class SquareMatrix : virtual public Matrix {
    public:
        enum SquareMatrixMethod {
            GaussMethod,
            LUMethod
        };
    public:
        /* Getters */
        [[nodiscard]] Scalar determinate(SquareMatrixMethod method);
    protected:
        /* Friends */
        friend class LinearEquations;
    };
    /////////////////////////////////////ColumnSquareMatrix//////////////////////////////////////////
    class ColumnSquareMatrix : public ColumnMatrix, public SquareMatrix {
    public:
        ColumnSquareMatrix();
        explicit ColumnSquareMatrix(size_t capacity);
        explicit ColumnSquareMatrix(const CStyleArray<Vector, Dynamic>& array);
        explicit ColumnSquareMatrix(CStyleArray<Vector, Dynamic>&& array) noexcept;
        ColumnSquareMatrix(const ColumnSquareMatrix& matrix);

        static ColumnSquareMatrix getUnitMatrix(size_t size);
    };
    /////////////////////////////////////RowSquareMatrix//////////////////////////////////////////
    class RowSquareMatrix : public RowMatrix, public SquareMatrix {
    public:
        RowSquareMatrix();
        explicit RowSquareMatrix(size_t capacity);
        explicit RowSquareMatrix(const CStyleArray<Vector, Dynamic>& array);
        explicit RowSquareMatrix(CStyleArray<Vector, Dynamic>&& array) noexcept;
        RowSquareMatrix(const RowSquareMatrix& matrix);

        static RowSquareMatrix getUnitMatrix(size_t size);
    };
}

#endif
