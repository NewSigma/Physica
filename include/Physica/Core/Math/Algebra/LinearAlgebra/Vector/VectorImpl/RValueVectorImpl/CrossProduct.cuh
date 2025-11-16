/*
 * Copyright 2023-2024 Weibo He.
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

#include "CrossProduct.h"

namespace Physica {
    template<Vector V1, Vector V2>
    class device_obj<CrossProduct<V1, V2>> : public device_obj<RValueVector<CrossProduct<V1, V2>>> {
        using host_obj = CrossProduct<V1, V2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    private:
        const device_obj<V1>& v1;
        const device_obj<V2>& v2;
    public:
        device_obj() = delete;
        __device__ device_obj(const device_obj<V1>& v1_, const device_obj<V2>& v2_) : v1(v1_), v2(v2_) {}
        /* Operations */
        __device__ void assign(Vector auto&& v) const {
            Base::assert_assign(v);
            v[0] = v1[1] * v2[2] - v1[2] * v2[1];
            v[1] = v1[2] * v2[0] - v1[0] * v2[2];
            v[2] = v1[0] * v2[1] - v1[1] * v2[0];
        }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getLength() const noexcept { return 3; }
    };
}

namespace Physica {
    template<Vector V1, Vector V2>
    class Traits<device_obj<CrossProduct<V1, V2>>> : public Traits<CrossProduct<V1, V2>> {};
}
