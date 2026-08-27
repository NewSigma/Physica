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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica {
    /**
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:236
     */
    template<Matrix M, Vector V>
    void applyHouseholder(const Scalar auto factor, const V& householder, M&& mat) {
        constexpr size_t BufferSize = householder.getSizeAtCompile() == Dynamic ? Dynamic : (householder.getSizeAtCompile() + 1);
        using T = V::ScalarType;
        using BufferVec = DenseVector<T, BufferSize>;
        assert(householder.getLength() + 1 == mat.getRow());
        auto copy = BufferVec(householder.getLength() + 1);
        copy[0] = 1;
        copy.tail(1) = householder;

        using BufferMat = DenseMatrix<T, MatrixMajor::Row, 1, decltype(mat * copy)::getSizeAtCompile()>;
        ((copy * -factor) * BufferMat(copy.hermite() * mat)).assign_add(mat);
    }

    template<Matrix M, Vector V>
    void applyHouseholder(M&& mat, const Scalar auto factor, const V& householder) {
        constexpr size_t BufferSize = householder.getSizeAtCompile() == Dynamic ? Dynamic : (householder.getSizeAtCompile() + 1);
        using T = V::ScalarType;
        using BufferVec = DenseVector<T, BufferSize>;
        assert(householder.getLength() + 1 == mat.getCol());
        auto copy = BufferVec(householder.getLength() + 1);
        copy[0] = 1;
        copy.tail(1) = householder;

        using BufferMat = DenseVector<T, decltype(mat * copy)::getSizeAtCompile()>;
        (BufferMat(mat * copy) * (copy.hermite() * -factor)).assign_add(mat);
    }

    template<Matrix M, Vector V>
    void applyHouseholder(const V& householder, M&& mat) {
        using T = V::ScalarType;
        if (householder.getLength() == 1) [[unlikely]]
            mat *= T(1) - householder.calc(0);
        else
            applyHouseholder(householder.calc(0), householder.tail(1), mat);
    }

    template<Matrix M, Vector V>
    void applyHouseholder(M&& mat, const V& householder) {
        using T = V::ScalarType;
        if (householder.getLength() == 1) [[unlikely]]
            mat *= T(1) - householder.calc(0);
        else
            applyHouseholder(mat, householder.calc(0), householder.tail(1));
    }
}
