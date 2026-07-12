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

#include "../CompactMatrix.cuh"
#include "Flatten.h"

namespace Physica {
    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    class device_obj<Flatten<M>> : public device_obj<CompactVector<Flatten<M>>> {
        using host_obj = Flatten<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactVector<host_obj>>;

        device_obj<M>& mat;
    public:
        device_obj(device_obj<CompactMatrix<M>>& mat_) : mat(mat_.getDerived()) {}
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

    template<Matrix M> requires(std::remove_cvref_t<M>::isCompact())
    __host__ __device__ auto device_obj<Flatten<M>>::data(this auto&& self) noexcept {
        return self.mat.data();
    }
}
