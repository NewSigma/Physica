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
                    result *= (*this)[i][i];
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

    ColumnSquareMatrix::ColumnSquareMatrix(size_t capacity)
            : Matrix(CStyleArray<Vector, Dynamic>(capacity), Column) {}

    ColumnSquareMatrix::ColumnSquareMatrix(const CStyleArray<Vector, Dynamic>& array)
            : Matrix(array, Column) {}

    ColumnSquareMatrix::ColumnSquareMatrix(CStyleArray<Vector, Dynamic>&& array) noexcept
            : Matrix(std::move(array), Column) {}

    ColumnSquareMatrix::ColumnSquareMatrix(const ColumnSquareMatrix& matrix) //NOLINT
            : Matrix(matrix), ColumnMatrix(matrix) {}

    ColumnSquareMatrix ColumnSquareMatrix::getUnitMatrix(size_t size) {
        const Numerical& _0 = BasicConst::getInstance().get_0();
        const Numerical& _1 = BasicConst::getInstance().get_1();
        CStyleArray<Vector, Dynamic> array(size, size);
        for(size_t i = 0; i < size; ++i) {
            CStyleArray<Numerical, Dynamic, Dynamic> array1(size, size);
            for(size_t j = 0; j < i; ++j)
                array1.allocate(Numerical(_0), j);
            array1.allocate(Numerical(_1), i);
            for(size_t j = i + 1; j < size; ++j)
                array1.allocate(Numerical(_0), j);
            array.allocate(Vector(std::move(array1)), i);
        }
        return ColumnSquareMatrix(std::move(array));
    }
    /////////////////////////////////////RowSquareMatrix//////////////////////////////////////////
    RowSquareMatrix::RowSquareMatrix() : Matrix(Row) {}

    RowSquareMatrix::RowSquareMatrix(size_t capacity)
            : Matrix(CStyleArray<Vector, Dynamic>(capacity), Row) {}

    RowSquareMatrix::RowSquareMatrix(const CStyleArray<Vector, Dynamic>& array)
            : Matrix(array, Row) {}

    RowSquareMatrix::RowSquareMatrix(CStyleArray<Vector, Dynamic>&& array) noexcept
            : Matrix(std::move(array), Row) {}

    RowSquareMatrix::RowSquareMatrix(const RowSquareMatrix& matrix) //NOLINT
            : Matrix(matrix), RowMatrix(matrix) {}

    RowSquareMatrix RowSquareMatrix::getUnitMatrix(size_t size) {
        const Numerical& _0 = BasicConst::getInstance().get_0();
        const Numerical& _1 = BasicConst::getInstance().get_1();
        CStyleArray<Vector, Dynamic> array(size, size);
        for(size_t i = 0; i < size; ++i) {
            CStyleArray<Numerical, Dynamic, Dynamic> array1(size, size);
            for(size_t j = 0; j < i; ++j)
                array1.allocate(Numerical(_0), j);
            array1.allocate(Numerical(_1), i);
            for(size_t j = i + 1; j < size; ++j)
                array1.allocate(Numerical(_0), j);
            array.allocate(Vector(std::move(array1)), i);
        }
        return RowSquareMatrix(std::move(array));
    }
}