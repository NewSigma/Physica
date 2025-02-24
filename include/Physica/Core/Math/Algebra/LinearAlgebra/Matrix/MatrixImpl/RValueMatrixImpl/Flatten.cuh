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

#include "Flatten.h"

namespace Physica {
    template<Matrix T>
    class device_obj<FlattenR<T>> : public device_obj<RValueVector<FlattenR<T>>> {
        using This = device_obj<FlattenR<T>>;

        const device_obj<T>& mat;
    public:
        using host_obj = FlattenR<T>;
        using Base = device_obj<RValueVector<host_obj>>;
        using typename Base::ScalarType;
    public:
        __host__ __device__ device_obj(const device_obj<T>& mat_) : mat(mat_) {}
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
    };

    template<Matrix T>
    __device__ auto device_obj<FlattenR<T>>::calc(size_t index) const -> ScalarType {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        return mat.calcFromMajorMinor(major, minor);
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<device_obj<FlattenR<T>>> : public Traits<FlattenR<T>> {};
}
