/*
 * Copyright 2020-2024 Weibo He.
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
            std::swap(biasOrder[i], biasOrder[matrix.partialPivoting(i)]);
            decompositionColumn(i);
        }
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    void PLUDecomposition<T, type, maxRow, maxColumn>::swap(PLUDecomposition& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        matrix.swap(obj.matrix);
        biasOrder.swap(biasOrder);
    }
    /*!
     * Apply LU Decomposition on a column of Matrix \from, save the result to Matrix \to.
     *
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:32
     */
    template<class T, int type, size_t maxRow, size_t maxColumn>
    void PLUDecomposition<T, type, maxRow, maxColumn>::decompositionColumn(size_t col) {
        const auto startAlphaIndex = col + 1;
        for (size_t j = 1; j < startAlphaIndex; ++j) { //Start from 1, unnecessary to handle j = 0
            T temp(matrix(j, col));
            for (size_t k = 0; k < j; ++k)
                temp -= matrix(j, k) * matrix(k, col);
            matrix(j, col) = std::move(temp);
        }

        const auto r = matrix.getRow();
        for (size_t j = startAlphaIndex; j < r; ++j) {
            T temp(matrix(j, col));
            for (size_t k = 0; k < col; ++k)
                temp -= matrix(j, k) * matrix(k, col);
            matrix(j, col) = temp / matrix(col, col);
        }
    }
}
