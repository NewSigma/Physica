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
    void applyHouseholder(const V& householder, M& mat) {
        using ScalarType = M::ScalarType;
        using BufferType = DenseVector<ScalarType, V::SizeAtCompile>;
        BufferType copy = householder;
        auto temp = ScalarType(1);
        temp.swap(copy[0]);
        const BufferType temp1 = copy * temp;
        mat -= temp1 * (copy.hermite() * mat).compute();
    }

    template<Matrix M, Vector V>
    void applyHouseholder(M& mat, const V& householder) {
        using ScalarType = M::ScalarType;
        using BufferType = DenseVector<ScalarType, V::SizeAtCompile>;
        BufferType copy = householder;
        ScalarType temp = ScalarType(1);
        temp.swap(copy[0]);
        using ProductType = decltype(mat * copy);
        using BufferType1 = DenseVector<ScalarType, ProductType::SizeAtCompile>;
        mat -= BufferType1(mat * copy) * (copy.hermite() * temp);
    }
}
