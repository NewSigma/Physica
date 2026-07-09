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
#include "Physica/Core/Parallel/ThreadBlock.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<Flatten<M>> : public device_obj<RValueVector<Flatten<M>>> {
        using host_obj = Flatten<M>;
        using This = device_obj<Flatten<M>>;
        using Base = device_obj<RValueVector<host_obj>>;

        const device_obj<M>& mat;
    protected:
        using typename Base::T;
    public:
        __host__ __device__ device_obj(const device_obj<M>& mat_) : mat(mat_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
    };

    template<Matrix M>
    __device__ auto device_obj<Flatten<M>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        return mat.calcFromMajorMinor(major, minor, block);
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<Flatten<M>>> : public Traits<Flatten<M>> {};
}
