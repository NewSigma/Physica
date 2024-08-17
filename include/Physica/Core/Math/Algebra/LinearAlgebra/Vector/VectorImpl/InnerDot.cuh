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

namespace Physica::Core {
    template<class VectorType1, class VectorType2>
    __device__ typename Internal::InnerDotImpl<VectorType1, VectorType2>::ResultType
    operator*(const device_obj<RValueVector<VectorType1>>& v1, const device_obj<RValueVector<VectorType2>>& v2) {
        using ResultType = typename Internal::InnerDotImpl<VectorType1, VectorType2>::ResultType;

        assert(v1.getLength() == v2.getLength() && "[Error]: Dimensions do not match");
        ResultType result = 0;
        for (size_t i = 0; i < v1.getLength(); ++i)
            result += v1.calc(i) * v2.calc(i);
        return result;
    }
}
