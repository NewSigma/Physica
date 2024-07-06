/*
 * Copyright 2021-2024 WeiBo He.
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

#include <cstdlib>
#include <utility>
#include <Physica/Macro.h>

namespace Physica::Core {
    template<class Matrix>
    class LUDecomposition {
        const Matrix& matrix;
    public:
        explicit LUDecomposition(const Matrix& m) : matrix(m) { assert(m.getRow() == m.getColumn()); }
        /* Operations */
        template<class MatrixOut>
        void decompositionColumn(MatrixOut& out, size_t column);
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return matrix.getRow(); }
        [[nodiscard]] size_t getColumn() const noexcept { return matrix.getRow(); }
        [[nodiscard]] const Matrix& getMatrix() const noexcept { return matrix; }
    };
    /*!
     * Apply LU Decomposition on a column of Matrix \from, save the result to Matrix \to.
     *
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:32
     */
    template<class Matrix>
    template<class MatrixOut>
    void LUDecomposition<Matrix>::decompositionColumn(MatrixOut& out, size_t column) {
        using ScalarType = typename Traits<Matrix>::ScalarType;
        const auto startAlphaIndex = column + 1;
        for (size_t j = 1; j < startAlphaIndex; ++j) {
            ScalarType temp(out(j, column));
            for (size_t k = 0; k < j; ++k)
                temp -= out(j, k) * out(k, column);
            out(j, column) = std::move(temp);
        }

        const auto r = out.getRow();
        for (size_t j = startAlphaIndex; j < r; ++j) {
            ScalarType temp(out(j, column));
            for (size_t k = 0; k < column; ++k)
                temp -= out(j, k) * out(k, column);
            out(j, column) = temp / out(column, column);
        }
    }
}
