/*
 * Copyright 2021 WeiBo He.
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

namespace Physica::Core {
    namespace Internal {
        template<class T>
        class Traits;
    }

    template<class Matrix>
    class LUDecomposition {
        const Matrix& matrix;
    public:
        explicit LUDecomposition(const Matrix& m) : matrix(m) { assert(m.getRow() == m.getColumn()); }
        /* Operations */
        template<class MatrixOut>
        void decompositionColumn(MatrixOut& out, size_t column);
        /* Getters */
        [[nodiscard]] const Matrix& getMatrix() const noexcept { return matrix; }
        [[nodiscard]] size_t getOrder() const noexcept { return matrix.getOrder(); }
    };
    /*!
     * Apply LU Decomposition on a column of Matrix \from, save the result to Matrix \to.
     *
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009.32
     */
    template<class Matrix>
    template<class MatrixOut>
    void LUDecomposition<Matrix>::decompositionColumn(MatrixOut& out, size_t column) {
        using ScalarType = typename Internal::Traits<Matrix>::ScalarType;
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
