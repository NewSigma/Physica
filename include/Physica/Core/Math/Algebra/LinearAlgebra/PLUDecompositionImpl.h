/*
 * Copyright 2020-2021 WeiBo He.
 *
 * This file is part of Physica.
 *
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
#ifndef PHYSICA_PLUDecompositionIMPL_H
#define PHYSICA_PLUDecompositionIMPL_H

#include "Matrix/MatrixOperation.h"

namespace Physica::Core {
    template<class T, int type, size_t maxRow, size_t maxColumn>
    PLUDecomposition<T, type, maxRow, maxColumn>::PLUDecomposition(const DenseMatrix<T, type, maxRow, maxColumn>& m)
            : matrix(m), biasOrder(reinterpret_cast<size_t*>(malloc(m.getRow() * sizeof(size_t)))) {
        const auto rank = matrix.getRow();
        for(size_t i = 0; i < rank; ++i)
            biasOrder[i] = i;
        for(size_t i = 0; i < rank; ++i) {
            std::swap(biasOrder[i], biasOrder[MatrixOperation<T, type, maxRow, maxColumn>::partialPivoting(matrix, i)]);
            decompositionColumn(matrix, i);
        }
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    PLUDecomposition<T, type, maxRow, maxColumn>::PLUDecomposition(DenseMatrix<T, type, maxRow, maxColumn>&& m) noexcept
            : matrix(std::move(m)), biasOrder(reinterpret_cast<size_t*>(malloc(m.getRow() * sizeof(size_t)))) {}

    template<class T, int type, size_t maxRow, size_t maxColumn>
    PLUDecomposition<T, type, maxRow, maxColumn>::PLUDecomposition(const PLUDecomposition& l)
            : matrix(l.matrix), biasOrder(reinterpret_cast<size_t*>(malloc(l.matrix.getRow() * sizeof(size_t)))) {
        Q_UNUSED(type)
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)
        const auto rank = l.matrix.getRow();
        for(size_t i = 0; i < rank; ++i) {
            biasOrder[i] = l.biasOrder[i];
        }
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    PLUDecomposition<T, type, maxRow, maxColumn>::PLUDecomposition(PLUDecomposition&& l) noexcept
            : matrix(std::move(l.matrix)), biasOrder(l.biasOrder) {
        Q_UNUSED(type)
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)
        l.biasOrder = nullptr;
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    PLUDecomposition<T, type, maxRow, maxColumn>& PLUDecomposition<T, type, maxRow, maxColumn>::operator= (
            const PLUDecomposition& l) {
        Q_UNUSED(type)
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)
        if(this != &l) {
            matrix = l.matrix;
            realloc(biasOrder, matrix.getRow() * sizeof(size_t));
        }
        return *this;
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    PLUDecomposition<T, type, maxRow, maxColumn>& PLUDecomposition<T, type, maxRow, maxColumn>::operator= (
            PLUDecomposition&& l) noexcept {
        Q_UNUSED(type)
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)
        this->~PLUDecomposition();
        matrix = std::move(l.matrix);
        biasOrder = l.biasOrder;
        l.biasOrder = nullptr;
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    PLUDecomposition<T, type, maxRow, maxColumn>::~PLUDecomposition() {
        Q_UNUSED(type)
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)
        free(biasOrder);
    }
    /*!
     * Apply LU Decomposition on a column of Matrix \from, save the result to Matrix \to.
     *
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.32
     */
    template<class T, int type, size_t maxRow, size_t maxColumn>
    void PLUDecomposition<T, type, maxRow, maxColumn>::decompositionColumn(size_t column) {
        Q_UNUSED(type)
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)
        const auto startAlphaIndex = column + 1;
        for (size_t j = 1; j < startAlphaIndex; ++j) {
            T temp(matrix(j, column));
            for (size_t k = 0; k < j; ++k)
                temp -= matrix(j, k) * matrix(k, column);
            matrix(j, column) = std::move(temp);
        }

        const auto r = matrix.getRow();
        for (size_t j = startAlphaIndex; j < r; ++j) {
            T temp(matrix(j, column));
            for (size_t k = 0; k < column; ++k)
                temp -= matrix(j, k) * matrix(k, column);
            matrix(j, column) = temp / matrix(column, column);
        }
    }
}

#endif