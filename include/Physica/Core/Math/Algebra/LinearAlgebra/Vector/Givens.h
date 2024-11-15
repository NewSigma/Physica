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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h>

namespace Physica::Core {
    /**
     * Construct a givens transformation that eleminate the element in \param vector at index \param j
     */
    template<class VectorType>
    auto givens(const LValueVector<VectorType>& vector, size_t i, size_t j) {
        using ScalarType = typename VectorType::ScalarType;
        using ResultType = Vector2D<ScalarType>;
        if constexpr (ScalarType::isComplex) {
            using RealType = typename ScalarType::RealType;
            using RealResultType = Vector2D<RealType>;

            ScalarType x_i = vector[i];
            ScalarType x_j = vector[j];
            const auto alpha = givens(RealResultType{x_i.real(), x_i.imag()}, 0, 1);
            const auto beta = givens(RealResultType{x_j.real(), x_j.imag()}, 0, 1);
            const auto theta = givens(RealResultType{x_i.norm(), x_j.norm()}, 0, 1);
            const auto factor = ScalarType(alpha[0] * beta[0] + alpha[1] * beta[1], alpha[1] * beta[0] - alpha[0] * beta[1]);
            return ResultType{theta[0], theta[1] * factor};
        }
        else {
            ScalarType x_i = vector[i];
            ScalarType x_j = vector[j];
            ScalarType rep_norm = reciprocal(sqrt(square(x_i) + square(x_j)));
            ScalarType cos = x_i * rep_norm;
            ScalarType sin = x_j * rep_norm;
            return ResultType{cos, sin};
        }
    }
    /**
     * Apply givens on left
     */
    template<class MatrixType>
    void applyGivens(const Vector2D<typename MatrixType::ScalarType>& givens, LValueMatrix<MatrixType>& mat, size_t i, size_t j) {
        using ScalarType = typename MatrixType::ScalarType;
        auto row_i = mat.row(i);
        auto row_j = mat.row(j);
        const size_t length = row_i.getLength();
        const ScalarType g0 = givens[0];
        const ScalarType g1 = givens[1];
        for (size_t k = 0; k < length; ++k) {
            const ScalarType temp1 = row_i[k];
            const ScalarType temp2 = row_j[k];
            row_i[k] = temp1 * g0 + temp2 * g1;
            row_j[k] = temp2 * g0 - temp1 * g1.conjugate();
        }
    }
    /**
     * Apply givens on right
     */
    template<class MatrixType>
    void applyGivens(LValueMatrix<MatrixType>& mat, const Vector2D<typename MatrixType::ScalarType>& givens, size_t i, size_t j) {
        using ScalarType = typename MatrixType::ScalarType;
        auto col_i = mat.col(i);
        auto col_j = mat.col(j);
        const size_t length = col_i.getLength();
        const ScalarType g0 = givens[0];
        const ScalarType g1 = givens[1];
        for (size_t k = 0; k < length; ++k) {
            const ScalarType temp1 = col_i[k];
            const ScalarType temp2 = col_j[k];
            col_i[k] = temp1 * g0 - temp2 * g1.conjugate();
            col_j[k] = temp1 * g1 + temp2 * g0;
        }
    }
}
