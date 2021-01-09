/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
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
    template<class T, DenseMatrixType type, size_t maxSize>
    SquareMatrix<T, type, maxSize>::SquareMatrix() : DenseMatrix<T, type, maxSize, maxSize>() {
        Q_UNUSED(type)
        Q_UNUSED(maxSize)
    }

    template<class T, DenseMatrixType type, size_t maxSize>
    SquareMatrix<T, type, maxSize>::SquareMatrix(size_t length) : DenseMatrix<T, type, maxSize, maxSize>(length) {
        Q_UNUSED(type)
        Q_UNUSED(maxSize)
    }

    template<class T, DenseMatrixType type, size_t maxSize>
    SquareMatrix<T, type, maxSize>::SquareMatrix(SquareMatrix&& matrix) noexcept
            : DenseMatrix<T, type, maxSize, maxSize>(std::move(matrix)) {
        Q_UNUSED(type)
        Q_UNUSED(maxSize)
    }

    template<class T, DenseMatrixType type, size_t maxSize>
    SquareMatrix<T, type, maxSize>& SquareMatrix<T, type, maxSize>::operator=(
            SquareMatrix<T, type, maxSize>&& m) noexcept {
        DenseMatrix<T, type, maxSize, maxSize>::operator=(std::move(m));
    }
    /*!
     * Note: This function will broke the origin matrix.
     *
     * Reference: MultiScalar Recipes in C++
     */
    template<class T, DenseMatrixType type, size_t maxSize>
    T SquareMatrix<T, type, maxSize>::determinate(DeterminateMethod method) {
        typedef MatrixOperation<T, type, maxSize, maxSize> Operation;

        const auto rank = Base::getRow();
        switch(rank) {
            case 1:
                return (*this)[0][0];
            case 2:
                return (*this)[0][0] * (*this)[1][1] - (*this)[1][0] * (*this)[0][1];
            default:
                MultiScalar result(BasicConst::getInstance()._1);
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
    }
    
    template<class T, DenseMatrixType type, size_t maxSize>
    SquareMatrix<T, type, maxSize> SquareMatrix<T, type, maxSize>::getUnitMatrix(size_t length) {
        Q_UNUSED(type)
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