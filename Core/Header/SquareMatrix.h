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
        /* Getters */
        [[nodiscard]] Numerical determinate(SquareMatrixMethod method);
    protected:
        inline void LUDecompositionColumn(size_t column);
        inline static Vector* unitMatrixVector(size_t n);
        /* Friends */
        friend class LinearEquations;
    };
    /////////////////////////////////////ColumnSquareMatrix//////////////////////////////////////////
    class ColumnSquareMatrix : public ColumnMatrix, public SquareMatrix {
    public:
        ColumnSquareMatrix();
        ColumnSquareMatrix(size_t column, size_t row);
        ColumnSquareMatrix(Vector* vectors, size_t length);

        static ColumnSquareMatrix unitColumnMatrix(size_t n) { return ColumnSquareMatrix(unitMatrixVector(n), n); }
    };
    /////////////////////////////////////RowSquareMatrix//////////////////////////////////////////
    class RowSquareMatrix : public RowMatrix, public SquareMatrix {
    public:
        RowSquareMatrix();
        RowSquareMatrix(size_t column, size_t row);
        RowSquareMatrix(Vector* vectors, size_t length);

        static RowSquareMatrix unitMatrix(size_t n) { return RowSquareMatrix(unitMatrixVector(n), n); }
    };
    /* Inline Implementations */
    inline void SquareMatrix::LUDecompositionColumn(size_t column) {
        const auto startAlphaIndex = column + 1;
        const auto rank = row();
        auto& matrix = *this;
        for(size_t j = 0; j < startAlphaIndex; ++j) {
            Numerical temp(matrix(j, column));
            for(size_t k = 0; k < j - 1; ++k)
                temp -= matrix(j, k) * matrix(k, column);
            matrix(j, column) = temp;
        }

        for(size_t j = startAlphaIndex; j < rank; ++j) {
            Numerical temp(matrix(j, column));
            for(size_t k = 0; k < j - 1; ++k)
                temp -= matrix(j, k) * matrix(k, column);
            matrix(j, column) = temp / matrix(column, column);
        }
    }
    //Implementation of unitColumnMatrix and unitRowMatrix. Should not be used elsewhere.
    inline Vector* SquareMatrix::unitMatrixVector(size_t n) {
        auto vectors = new Vector[n];
        auto& _0 = BasicConst::getInstance().get_0();
        auto& _1 = BasicConst::getInstance().get_1();
        for(size_t i = 0; i < n; ++i) {
            auto& vector = vectors[i];
            vector.initVector(n);
            for(size_t j = 0; j < i; ++j)
                vector[j] = _0;
            vector[i] = _1;
            for(size_t j = i + 1; j < n; ++j)
                vector[j] = _0;
        }
        return vectors;
    }
}

#endif
