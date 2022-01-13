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

#include "LValueVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/LValueMatrix.h"

namespace Physica::Core {
    namespace Internal {
        template<class VectorType, bool isComplex>
        struct GivensImpl {
            using ScalarType = typename VectorType::ScalarType;
            using ResultType = Vector<ScalarType, 2, 2>;

            static ResultType run(const LValueVector<VectorType>& vector, size_t i, size_t j) {
                ScalarType x_i = vector[i];
                ScalarType x_j = vector[j];
                ScalarType rep_norm = reciprocal(sqrt(square(x_i) + square(x_j)));
                ScalarType cos = x_i * rep_norm;
                ScalarType sin = x_j * rep_norm;
                return ResultType{cos, sin};
            }
        };

        template<class VectorType>
        struct GivensImpl<VectorType, true> {
            using ScalarType = typename VectorType::ScalarType;
            using ResultType = Vector<ScalarType, 2, 2>;
            using RealType = typename ScalarType::RealType;
            using RealResultType = Vector<RealType, 2, 2>;
            using RealGivens = GivensImpl<RealResultType, false>;

            static ResultType run(const LValueVector<VectorType>& vector, size_t i, size_t j) {
                ScalarType x_i = vector[i];
                ScalarType x_j = vector[j];
                const auto alpha = RealGivens::run(RealResultType{x_i.getReal(), x_i.getImag()}, 0, 1);
                const auto beta = RealGivens::run(RealResultType{x_j.getReal(), x_j.getImag()}, 0, 1);
                const auto theta = RealGivens::run(RealResultType{x_i.norm(), x_j.norm()}, 0, 1);
                const ScalarType factor = ScalarType(alpha[0] * beta[0] + alpha[1] * beta[1], alpha[1] * beta[0] - alpha[0] * beta[1]);
                return ResultType{theta[0], theta[1] * factor};
            }
        };
    }
    /**
     * Construct a givens transformation that eleminate the element in \param vector at index \param j
     */
    template<class VectorType>
    Vector<typename VectorType::ScalarType, 2, 2> givens(const LValueVector<VectorType>& vector, size_t i, size_t j) {
        using ScalarType = typename VectorType::ScalarType;
        return Internal::GivensImpl<VectorType, ScalarType::isComplex>::run(vector, i, j);
    }
    /**
     * Apply givens on left
     */
    template<class MatrixType>
    void applyGivens(const Vector<typename MatrixType::ScalarType, 2>& givens, LValueMatrix<MatrixType>& mat, size_t i, size_t j) {
        auto row_i = mat.row(i);
        auto row_j = mat.row(j);
        const size_t length = row_i.getLength();
        for (size_t k = 0; k < length; ++k) {
            auto temp1 = row_i[k];
            auto temp2 = row_j[k];
            row_i[k] = temp1 * givens[0] + temp2 * givens[1];
            row_j[k] = -temp1 * givens[1].conjugate() + temp2 * givens[0];
        }
    }
    /**
     * Apply givens on right
     */
    template<class MatrixType>
    void applyGivens(LValueMatrix<MatrixType>& mat, const Vector<typename MatrixType::ScalarType, 2>& givens, size_t i, size_t j) {
        auto col_i = mat.col(i);
        auto col_j = mat.col(j);
        const size_t length = col_i.getLength();
        for (size_t k = 0; k < length; ++k) {
            auto temp1 = col_i[k];
            auto temp2 = col_j[k];
            col_i[k] = temp1 * givens[0] - temp2 * givens[1].conjugate();
            col_j[k] = temp1 * givens[1] + temp2 * givens[0];
        }
    }
}
