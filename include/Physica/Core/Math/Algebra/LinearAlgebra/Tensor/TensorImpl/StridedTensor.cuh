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
#include "StridedTensor.h"

namespace Physica {
    template<class Derived>
    class device_obj<StridedTensor<Derived>> : public device_obj<LValueTensor<Derived>> {
        using host_obj = StridedTensor<Derived>;
        using Base = device_obj<LValueTensor<Derived>>;
        using This = device_obj<host_obj>;
    public:
        using typename Base::IndexType;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto getStrides() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getStride(size_t dim) const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, const IndexType& index) noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, std::integral auto... dims) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static IndexType getStrideAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getStrideAtCompile(size_t dim) noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "StridedTensorImpl/StridedTensorImpl.cuh"
