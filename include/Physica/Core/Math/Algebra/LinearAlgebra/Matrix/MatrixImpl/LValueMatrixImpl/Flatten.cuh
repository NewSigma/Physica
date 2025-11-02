/*
 * Copyright 2023-2025 Weibo He.
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
    class device_obj<FlattenL<T>> : public device_obj<LValueVector<FlattenL<T>>> {
        using host_obj = FlattenL<T>;
        using This = device_obj<host_obj>;

        const device_obj<T>& mat;
    public:
        using Base = device_obj<LValueVector<host_obj>>;
        using typename Base::ScalarType;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    public:
        __host__ __device__ device_obj(const device_obj<LValueMatrix<T>>& mat_) : mat(mat_) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const;
    };

    template<Matrix T>
    __host__ __device__ auto device_obj<FlattenL<T>>::data_ptr(size_t index) -> PtrTy {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        const size_t row = MatrixOption::rowFromMajorMinor<T>(major, minor);
        const size_t col = MatrixOption::colFromMajorMinor<T>(major, minor);
        return mat.data_ptr(row, col);
    }

    template<Matrix T>
    __host__ __device__ auto device_obj<FlattenL<T>>::data_ptr(size_t index) const -> ConstPtrTy {
        return const_cast<This&>(*this).data_ptr(index);
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<device_obj<FlattenL<T>>> : public Traits<FlattenL<T>> {};
}
