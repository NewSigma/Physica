/*
 * Copyright 2023-2026 Weibo He.
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

#include "../ContinuousMatrix.cuh"
#include "Flatten.h"

namespace Physica {
    template<Matrix M>
    class device_obj<FlattenC<M>> : public device_obj<ContinuousVector<FlattenC<M>>> {
        using host_obj = FlattenC<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<ContinuousVector<host_obj>>;

        device_obj<M>& mat;
    public:
        device_obj(device_obj<ContinuousMatrix<M>>& mat_) : mat(mat_.getDerived()) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
        [[nodiscard]] __host__ __device__ auto data(this auto&&) noexcept;
    };

    template<Matrix M>
    __host__ __device__ auto device_obj<FlattenC<M>>::data(this auto&& self) noexcept {
        return self.mat.data();
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<FlattenC<M>>> : public Traits<FlattenC<M>> {};
}
