/*
 * Copyright 2024 Weibo He.
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

#include "../RValueVector.cuh"

namespace Physica {
    template<Vector T1, Vector T2>
    __device__ auto operator*(const T1& v1, const T2& v2) requires(CUDA<T1> && CUDA<T2>) {
        using ScalarType = InnerDot<T1, T2>::ScalarType;

        assert(v1.getLength() == v2.getLength() && "[Error]: Dimensions do not match");
        ScalarType result = 0;
        for (size_t i = 0; i < v1.getLength(); ++i)
            result += v1.calc(i) * v2.calc(i);
        return result;
    }
}
