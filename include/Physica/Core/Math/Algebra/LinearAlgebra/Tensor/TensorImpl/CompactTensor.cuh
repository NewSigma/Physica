/*
 * Copyright 2026 Weibo He.
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

#include "LValueTensor.cuh"
#include "CompactTensor.h"

namespace Physica {
    template<class Derived>
    class device_obj<CompactTensor<Derived>> : public device_obj<LValueTensor<Derived>> {
        using host_obj = CompactTensor<Derived>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueTensor<Derived>>;
    public:
        using typename Base::IndexType;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto data(this auto&& self) noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, const IndexType& index) noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, std::integral auto... dims) noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "CompactTensorImpl/CompactTensorImpl.cuh"
