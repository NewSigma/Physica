/*
 * Copyright 2020-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"

namespace Physica {
    /**
     * Construct a givens transformation that eleminate the element in \param vector at index \param j
     *
     * Reference:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:240-244
     */
    template<Vector V>
    Vector2D<typename V::ScalarType> givens(const V& vector, size_t i, size_t j) {
        static_assert(!ReverseDiff<V>, "[Error]: Not implemented, wait for P2014Rx");
        using T = V::ScalarType;
        using Tr = T::RealType;
        using ResultType = Vector2D<T>;
        using Vector2Dr = Vector2D<Tr>;
        T x_i = vector.calc(i);
        T x_j = vector.calc(j);
        ResultType result;
        if constexpr (T::isComplex) {
            const auto alpha = givens(Vector2Dr{x_i.real(), x_i.imag()}, 0, 1);
            const auto beta = givens(Vector2Dr{x_j.real(), x_j.imag()}, 0, 1);
            const auto theta = givens(Vector2Dr{x_i.norm(), x_j.norm()}, 0, 1);
            const auto factor = T(alpha[0] * beta[0] + alpha[1] * beta[1], alpha[1] * beta[0] - alpha[0] * beta[1]);
            const auto theta1 = theta[1] * factor;
            result = ResultType{theta[0], theta1};
        }
        else {
            if (x_j.isZero())
                return ResultType{Tr(x_i.isPositive() ? 1.0 : -1.0), 0};

            T rep_norm = reciprocal(sqrt(square(x_i) + square(x_j)));
            T cos = x_i * rep_norm;
            T sin = x_j * rep_norm;
            result = ResultType{cos, sin};
        }
        return result;
    }
    /**
     * Apply givens on left
     */
    template<Matrix M>
    void applyGivens(const Vector2D<typename M::ScalarType>& givens, LValueMatrix<M>& mat, size_t i, size_t j) {
        using T = M::ScalarType;
        auto row_i = mat.row(i);
        auto row_j = mat.row(j);
        const size_t length = row_i.getLength();
        const T g0 = givens[0];
        const T g1 = givens[1];
        for (size_t k = 0; k < length; ++k) {
            const T temp1 = row_i[k];
            const T temp2 = row_j[k];
            row_i[k] = temp1 * g0 + temp2 * g1;
            row_j[k] = temp2 * g0 - temp1 * g1.conjugate();
        }
    }
    /**
     * Apply givens on right
     */
    template<Matrix M>
    void applyGivens(LValueMatrix<M>& mat, const Vector2D<typename M::ScalarType>& givens, size_t i, size_t j) {
        using T = M::ScalarType;
        auto col_i = mat.col(i);
        auto col_j = mat.col(j);
        const size_t length = col_i.getLength();
        const T g0 = givens[0];
        const T g1 = givens[1];
        for (size_t k = 0; k < length; ++k) {
            const T temp1 = col_i[k];
            const T temp2 = col_j[k];
            col_i[k] = temp1 * g0 - temp2 * g1.conjugate();
            col_j[k] = temp1 * g1 + temp2 * g0;
        }
    }
}
