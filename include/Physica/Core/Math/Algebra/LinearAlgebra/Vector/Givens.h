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

namespace Physica::Core {
    template<class ScalarType>
    class Givens : public Vector<ScalarType, 2, 2> {};
    /**
     * Construct a givens transformation that eleminate the element in \param vector at index \param j
     */
    template<class VectorType>
    Givens<typename VectorType::ScalarType> givens(const VectorType& vector, size_t i, size_t j) {
        using ScalarType = typename VectorType::ScalarType;
        using ResultType = Givens<ScalarType>;
        ScalarType x_i = vector[i];
        ScalarType x_j = vector[j];
        ScalarType rep_norm = reciprocal(sqrt(square(x_i) + square(x_j)));
        ScalarType cos = x_i * rep_norm;
        ScalarType sin = x_j * rep_norm;
        return ResultType{cos, sin};
    }
    /**
     * Apply givens on left
     */
    template<class MatrixType, class T>
    void applyGivens(const Givens<T>& givens, MatrixType& mat, size_t i, size_t j) {
        using ScalarType = typename MatrixType::ScalarType;
        auto row_i = mat.row(i);
        auto row_j = mat.row(j);
        const size_t length = row_i.getLength();
        for (size_t k = 0; k < length; ++k) {
            auto temp1 = row_i[k];
            auto temp2 = row_j[k];
            row_i[k] = temp1 * givens[0] - temp2 * givens[1];
            row_j[k] = temp1 * givens[1] + temp2 * givens[0];
        }
    }
    /**
     * Apply givens on right
     */
    template<class MatrixType, class T>
    void applyGivens(MatrixType& mat, const Givens<T>& givens, size_t i, size_t j) {
        using ScalarType = typename MatrixType::ScalarType;
        auto col_i = mat.col(i);
        auto col_j = mat.col(j);
        const size_t length = col_i.getLength();
        for (size_t k = 0; k < length; ++k) {
            auto temp1 = col_i[k];
            auto temp2 = col_j[k];
            col_i[k] = temp1 * givens[0] + temp2 * givens[1];
            col_j[k] = -temp1 * givens[1] + temp2 * givens[0];
        }
    }
}
