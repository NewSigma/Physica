/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "SquareMatrix.h"

namespace Physica::Core {
    /////////////////////////////////////SquareMatrix//////////////////////////////////////////
    /*!
     * SquareMatrix is the matrix whose row() equals to column().
     * If the column of a SquareMatrix less than its row, out of bounder visiting will happen during
     * the latter calculation. If the column of a SquareMatrix more than its row,
     * the values of unnecessary columns will may also be changed.
     */
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
            case LUMethod:
                for(size_t i = 0; i < rank; ++i) {
                    partialPivoting(i);
                    LUDecompositionColumn(i);
                    result *= vectors[i][i];
                }
                break;
            default:;
        }
        return result;
    }
    /////////////////////////////////////ColumnSquareMatrix//////////////////////////////////////////
    ColumnSquareMatrix::ColumnSquareMatrix() : Matrix(Column) {}

    ColumnSquareMatrix::ColumnSquareMatrix(size_t column, size_t row)
            : Matrix(new Vector[column], column, Column), ColumnMatrix(column, row) {}

    ColumnSquareMatrix::ColumnSquareMatrix(Vector* vectors, size_t length) : Matrix(vectors, length, Column) {}
    /////////////////////////////////////RowSquareMatrix//////////////////////////////////////////
    RowSquareMatrix::RowSquareMatrix() : Matrix(Row) {}

    RowSquareMatrix::RowSquareMatrix(size_t column, size_t row)
            : Matrix(new Vector[row], row, Row), RowMatrix(column, row) {}

    RowSquareMatrix::RowSquareMatrix(Vector* vectors, size_t length) : Matrix(vectors, length, Row) {}
}