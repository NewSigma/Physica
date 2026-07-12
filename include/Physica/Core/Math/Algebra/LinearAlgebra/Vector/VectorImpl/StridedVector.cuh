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

#include "LValueVector.cuh"
#include "StridedVector.h"

namespace Physica {
    template<class Derived>
    class device_obj<StridedVector<Derived>> : public device_obj<LValueVector<Derived>> {
        using host_obj = StridedVector<Derived>;
        using Base = device_obj<LValueVector<Derived>>;
        using This = device_obj<host_obj>;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getStride() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t index) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getStrideAtCompile() noexcept { return host_obj::getStrideAtCompile(); }
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };

    template<class Derived>
    __host__ __device__ constexpr size_t device_obj<StridedVector<Derived>>::getStride() const noexcept {
        if constexpr (Derived::getStrideAtCompile() != Dynamic)
            return Derived::getStrideAtCompile();
        else
            return Base::getDerived().getStride();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedVector<Derived>>::data_handle() noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedVector<Derived>>::data_handle() const noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedVector<Derived>>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.getLength());
        return self.data_handle() + index * self.getStride();
    }
}
