/*
 * Copyright 2024-2026 Weibo He.
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
    template<Vector V1, Vector V2>
    __device__ auto operator*(const V1& v1, const V2& v2) requires(DeviceObj<V1> && DeviceObj<V2>) {
        using T1 = V1::ScalarType;
        using T2 = V2::ScalarType;
        using T = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;

        assert(v1.getLength() == v2.getLength() && "[Error]: Dimensions do not match");
        T result = 0;
        for (size_t i = 0; i < v1.getLength(); ++i) {
            if constexpr (std::same_as<T1, T2>)
                result = fma(v1.calc(i), v2.calc(i), result);
            else
                result += v1.calc(i) * v2.calc(i);
        }
        return result;
    }
}
