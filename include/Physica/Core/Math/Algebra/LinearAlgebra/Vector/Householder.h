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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    /**
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:236
     */
    template<Matrix M, Vector V, Scalar T>
    void applyHouseholder(const T factor, const V& householder, M&& mat) {
        constexpr size_t BufferSize = V::SizeAtCompile == Dynamic ? Dynamic : (V::SizeAtCompile + 1);
        using BufferType = DenseVector<T, BufferSize>;
        assert(householder.getLength() + 1 == mat.getRow());
        auto copy = BufferType(householder.getLength() + 1);
        copy[0] = 1;
        copy.tail(1) = householder;

        const BufferType temp1 = copy * factor;
        mat -= temp1 * (copy.hermite() * mat).compute();
    }

    template<Matrix M, Vector V, Scalar T>
    void applyHouseholder(M&& mat, const T factor, const V& householder) {
        constexpr size_t BufferSize = V::SizeAtCompile == Dynamic ? Dynamic : (V::SizeAtCompile + 1);
        using BufferType = DenseVector<T, BufferSize>;
        assert(householder.getLength() + 1 == mat.getCol());
        auto copy = BufferType(householder.getLength() + 1);
        copy[0] = 1;
        copy.tail(1) = householder;

        using BufferType1 = DenseVector<T, decltype(mat * copy)::SizeAtCompile>;
        mat -= BufferType1(mat * copy) * (copy.hermite() * factor);
    }

    template<Matrix M, Vector V>
    void applyHouseholder(const V& householder, M&& mat) {
        applyHouseholder(householder[0], householder.tail(1), mat);
    }

    template<Matrix M, Vector V>
    void applyHouseholder(M&& mat, const V& householder) {
        applyHouseholder(mat, householder[0], householder.tail(1));
    }
}
