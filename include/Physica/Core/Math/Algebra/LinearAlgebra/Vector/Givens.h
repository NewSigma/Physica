/*
 * Copyright 2020-2026 Weibo He.
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
    auto givens(const V& vector, size_t i, size_t j) noexcept {
        static_assert(!ReverseDiff<V>, "[Error]: Not implemented, wait for P2014Rx");
        using T = V::ScalarType;
        using Tr = T::RealType;
        using ResultType = Vector2D<T>;
        using Vector2Dr = Vector2D<Tr>;
        T x_i = vector.calc(i);
        T x_j = vector.calc(j);
        if constexpr (T::isComplex()) {
            const auto alpha = givens(Vector2Dr{x_i.real(), x_i.imag()}, 0, 1);
            const auto beta = givens(Vector2Dr{x_j.real(), x_j.imag()}, 0, 1);
            const auto theta = givens(Vector2Dr{x_i.norm(), x_j.norm()}, 0, 1);
            const auto factor = T(fma(alpha[0], beta[0], alpha[1] * beta[1]), fma(alpha[1], beta[0], -alpha[0] * beta[1]));
            const auto theta1 = theta[1] * factor;
            return ResultType{theta[0], theta1};
        }
        else {
            if (x_j.isZero())
                return ResultType{Tr(x_i.isPositive() ? 1.0 : -1.0), Tr(0)};

            T rep_norm = reciprocal(sqrt(fma(x_i, x_i, square(x_j))));
            T cos = x_i * rep_norm;
            T sin = x_j * rep_norm;
            return ResultType{cos, sin};
        }
    }
    /**
     * Givens * Matrix
     */
    template<Matrix M>
    void applyGivens(const Vector2D<typename std::remove_cvref_t<M>::ScalarType>& givens, M&& mat, size_t i, size_t j) noexcept {
        using T = std::remove_cvref_t<M>::ScalarType;
        const auto cosine = givens[0].real();
        const auto sine = givens[1];
        auto row_i = mat.row(i);
        auto row_j = mat.row(j);

        const size_t length = row_i.getLength();
        size_t k = 0;
        if constexpr (row_i.isCompact()) {
            using Pack = BestPacket<T, mat.getColAtCompile()>::Type;
            constexpr int Size = Pack::size();
            auto it = zip(row_i.view(), row_j.view()).begin();
            for (; k < length / Size * Size; k += Size) {
                auto [it_i, it_j] = it + k;
                auto pack1 = it_i.template load<Size>();
                auto pack2 = it_j.template load<Size>();
                it_i.store(fma(pack1, Pack(cosine), pack2 * Pack(sine)));
                it_j.store(fma(pack2, Pack(cosine), pack1 * Pack(-sine.conjugate())));
            }
        }

        for (; k < length; ++k) {
            const T temp1 = row_i[k];
            const T temp2 = row_j[k];
            row_i[k] = fma(temp1, cosine, temp2 * sine);
            row_j[k] = fma(temp2, cosine, temp1 * (-sine.conjugate()));
        }
    }
    /**
     * Matrix * Givens
     */
    template<Matrix M>
    void applyGivens(M&& mat, const Vector2D<typename std::remove_cvref_t<M>::ScalarType>& givens, size_t i, size_t j) noexcept {
        using T = std::remove_cvref_t<M>::ScalarType;
        const auto cosine = givens[0].real();
        const auto sine = givens[1];
        auto col_i = mat.col(i);
        auto col_j = mat.col(j);

        const size_t length = col_i.getLength();
        size_t k = 0;
        if constexpr (col_i.isCompact()) {
            using Pack = BestPacket<T, mat.getColAtCompile()>::Type;
            constexpr int Size = Pack::size();
            auto it = zip(col_i.view(), col_j.view()).begin();
            for (; k < length / Size * Size; k += Size) {
                auto [it_i, it_j] = it + k;
                auto pack1 = it_i.template load<Size>();
                auto pack2 = it_j.template load<Size>();
                it_i.store(fma(pack1, Pack(cosine), pack2 * Pack(-sine.conjugate())));
                it_j.store(fma(pack1, Pack(sine), pack2 * Pack(cosine)));
            }
        }

        for (; k < length; ++k) {
            const T temp1 = col_i[k];
            const T temp2 = col_j[k];
            col_i[k] = fma(temp1, cosine, temp2 * (-sine.conjugate()));
            col_j[k] = fma(temp1, sine, temp2 * cosine);
        }
    }
    /**
     * Givens * Vector
     */
    template<Vector V>
    void applyGivensCol(const Vector2D<typename std::remove_cvref_t<V>::ScalarType>& givens, V&& vec, size_t i, size_t j) noexcept {
        using T = std::remove_cvref_t<V>::ScalarType;
        const auto cosine = givens[0].real();
        const auto sine = givens[1];
        const T temp1 = vec[i];
        const T temp2 = vec[j];
        vec[i] = fma(temp1, cosine, temp2 * sine);
        vec[j] = fma(temp2, cosine, temp1 * (-sine.conjugate()));
    }
    /**
     * Vector * Givens
     */
    template<Vector V>
    void applyRowGivens(V&& vec, const Vector2D<typename std::remove_cvref_t<V>::ScalarType>& givens, size_t i, size_t j) noexcept {
        using T = std::remove_cvref_t<V>::ScalarType;
        const auto cosine = givens[0].real();
        const auto sine = givens[1];
        const T temp1 = vec[i];
        const T temp2 = vec[j];
        vec[i] = fma(temp1, cosine, temp2 * (-sine.conjugate()));
        vec[j] = fma(temp1, sine, temp2 * cosine);
    }
}
