/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include <Core/Header/LUDecomposition.h>
#include "Core/Header/SquareMatrix.h"

namespace Physica::Core {
    /////////////////////////////////////SquareMatrix//////////////////////////////////////////
    /*!
     * SquareMatrix is the matrix whose row() equals to column().
     * If the column of a SquareMatrix less than its row, out of bounder visiting will happen during
     * the latter calculation. If the column of a SquareMatrix more than its row,
     * the values of unnecessary columns will may also be changed.
     */
    void SquareMatrix::toUnit() {
        auto& _0 = BasicConst::getInstance().get_0();
        auto& _1 = BasicConst::getInstance().get_1();
        for(size_t i = 0; i < length; ++i) {
            for(size_t j = 0; j < i; ++j)
                vectors[i][j] = _0;
            vectors[i][i] = _1;
            for(size_t j = i + 1; j < length; ++j)
                vectors[i][j] = _0;
        }
    }
    /*!
     * Note: This function will broke the origin matrix.
     *
     * Reference: Numerical Recipes in C++
     */
    Numerical SquareMatrix::determinate(SquareMatrixMethod method) {
        const auto rank = row();
        Numerical result(BasicConst::getInstance().get_1());
        switch(method) {
            case GaussMethod:
                for(size_t i = 0; i < rank; ++i) {
                    upperEliminate(i);
                    lowerEliminate(i);
                    result *= vectors[i][i];
                }
                break;
            case LUMethod: {
                LUDecomposition lu(*this);
                const Matrix& m = lu.getMatrix();
                for(size_t i = 0; i < rank; ++i)
                    result *= m(i, i);
            }
                break;
            default:;
        }
        return result;
    }
    /////////////////////////////////////ColumnSquareMatrix//////////////////////////////////////////
    ColumnSquareMatrix::ColumnSquareMatrix() : Matrix(Column) {}

    ColumnSquareMatrix::ColumnSquareMatrix(size_t size)
            : Matrix(reinterpret_cast<Vector*>(malloc(size * sizeof(Vector))), size, Column)
            , ColumnMatrix(size, size) {}

    ColumnSquareMatrix::ColumnSquareMatrix(Vector* vectors, size_t length) : Matrix(vectors, length, Column) {}

    ColumnSquareMatrix::ColumnSquareMatrix(const ColumnSquareMatrix& matrix) //NOLINT
            : Matrix(matrix), ColumnMatrix(matrix) {}
    /////////////////////////////////////RowSquareMatrix//////////////////////////////////////////
    RowSquareMatrix::RowSquareMatrix() : Matrix(Row) {}

    RowSquareMatrix::RowSquareMatrix(size_t size)
            : Matrix(reinterpret_cast<Vector*>(malloc(size * sizeof(Vector))), size, Row)
            , RowMatrix(size, size) {}

    RowSquareMatrix::RowSquareMatrix(Vector* vectors, size_t length) : Matrix(vectors, length, Row) {}

    RowSquareMatrix::RowSquareMatrix(const RowSquareMatrix& matrix) //NOLINT
            : Matrix(matrix), RowMatrix(matrix) {}
}