/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SQUAREMATRIXIMPL_H
#define PHYSICA_SQUAREMATRIXIMPL_H

#include "Physica/Core/Math/Algebra/LinearAlgebra/LUDecomposition.h"
#include "MatrixOperation.h"

namespace Physica::Core {
    /*!
     * SquareMatrix is the matrix whose row() equals to column().
     * If the column of a SquareMatrix less than its row, out of bounder visiting will happen during
     * the latter calculation. If the column of a SquareMatrix more than its row,
     * the values of unnecessary columns will may also be changed.
     */
    template<class T, MatrixType type, size_t maxSize>
    SquareMatrix<T, type, maxSize>::SquareMatrix() : Matrix<T, type, maxSize, maxSize>() {}

    template<class T, MatrixType type, size_t maxSize>
    SquareMatrix<T, type, maxSize>::SquareMatrix(size_t length) : Matrix<T, type, maxSize, maxSize>(length) {}

    template<class T, MatrixType type, size_t maxSize>
    SquareMatrix<T, type, maxSize>::SquareMatrix(SquareMatrix&& matrix) noexcept
            : Matrix<T, type, maxSize, maxSize>(std::move(matrix)) {}

    template<class T, MatrixType type, size_t maxSize>
    SquareMatrix<T, type, maxSize>& SquareMatrix<T, type, maxSize>::operator=(
            SquareMatrix<T, type, maxSize>&& m) noexcept {
        Matrix<T, type, maxSize, maxSize>::operator=(std::move(m));
    }
    /*!
     * Note: This function will broke the origin matrix.
     *
     * Reference: MultiScalar Recipes in C++
     */
    template<class T, MatrixType type, size_t maxSize>
    T SquareMatrix<T, type, maxSize>::determinate(DeterminateMethod method) {
        typedef MatrixOperation<T, type, maxSize, maxSize> Operation;

        const auto rank = Matrix<T, type, maxSize, maxSize>::getRow();
        MultiScalar result(BasicConst::getInstance().get_1());
        switch(method) {
            case GaussMethod:
                for(size_t i = 0; i < rank; ++i) {
                    Operation::upperEliminate(*this, i);
                    Operation::lowerEliminate(*this, i);
                    result *= (*this)[i][i];
                }
                break;
            case LUMethod: {
                LUDecomposition<T, type, maxSize, maxSize> lu(*this);
                const auto& m = lu.getMatrix();
                for(size_t i = 0; i < rank; ++i)
                    result *= m(i, i);
            }
                break;
            default:;
        }
        return result;
    }
    
    template<class T, MatrixType type, size_t maxSize>
    SquareMatrix<T, type, maxSize> SquareMatrix<T, type, maxSize>::getUnitMatrix(size_t length) {
        SquareMatrix result(length);
        for(size_t i = 0; i < length; ++i) {
            Vector<T, maxSize> vector(length);
            for(size_t j = 0; j < i; ++j)
                vector.allocate(T::getZero(), j);
            vector.allocate(T::getOne(), i);
            for(size_t j = i + 1; j < length; ++j)
                vector.allocate(T::getZero(), j);
            result.allocate(std::move(vector), i);
        }
        return result;
    }
}

#endif