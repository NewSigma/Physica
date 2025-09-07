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
    template<Matrix M, Vector V>
    void applyHouseholder(const V& householder, M&& mat) {
        using T = std::remove_cvref<M>::type::ScalarType;
        using BufferType = DenseVector<T, V::SizeAtCompile>;
        assert(householder.getLength() == mat.getRow());
        BufferType copy = householder;

        const T factor = copy[0];
        copy[0] = T(1);
        const BufferType temp1 = copy * factor;
        mat -= temp1 * (copy.hermite() * mat).compute();
    }

    template<Matrix M, Vector V>
    void applyHouseholder(M&& mat, const V& householder) {
        using T = std::remove_cvref<M>::type::ScalarType;
        using BufferType = DenseVector<T, V::SizeAtCompile>;
        assert(householder.getLength() == mat.getCol());
        BufferType copy = householder;

        using ProductType = decltype(mat * copy);
        using BufferType1 = DenseVector<T, ProductType::SizeAtCompile>;
        const T factor = copy[0];
        copy[0] = T(1);
        mat -= BufferType1(mat * copy) * (copy.hermite() * factor);
    }
}
