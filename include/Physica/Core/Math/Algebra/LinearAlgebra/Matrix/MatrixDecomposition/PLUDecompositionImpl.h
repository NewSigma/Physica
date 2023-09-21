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
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixOperation.h"

namespace Physica::Core {
    template<class T, int type, size_t maxRow, size_t maxColumn>
    PLUDecomposition<T, type, maxRow, maxColumn>::PLUDecomposition(MatrixType m) {
        compute(std::move(m));
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    PLUDecomposition<T, type, maxRow, maxColumn>&
    PLUDecomposition<T, type, maxRow, maxColumn>::operator=(PLUDecomposition obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    void PLUDecomposition<T, type, maxRow, maxColumn>::compute(MatrixType m) {
        matrix = std::move(m);
        biasOrder.resize(matrix.getRow());
        const auto rank = matrix.getRow();
        for(size_t i = 0; i < rank; ++i)
            biasOrder[i] = i;
        for(size_t i = 0; i < rank; ++i) {
            std::swap(biasOrder[i], biasOrder[MatrixOperation<T, type, maxRow, maxColumn>::partialPivoting(matrix, i)]);
            decompositionColumn(i);
        }
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    void PLUDecomposition<T, type, maxRow, maxColumn>::swap(PLUDecomposition& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        matrix.swap(obj.matrix);
        biasOrder.swap(biasOrder);
    }
    /*!
     * Apply LU Decomposition on a column of Matrix \from, save the result to Matrix \to.
     *
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009:32
     */
    template<class T, int type, size_t maxRow, size_t maxColumn>
    void PLUDecomposition<T, type, maxRow, maxColumn>::decompositionColumn(size_t column) {
        const auto startAlphaIndex = column + 1;
        for (size_t j = 1; j < startAlphaIndex; ++j) { //Start from 1, unnecessary to handle j = 0
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
