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

#include "LValueMatrix.cuh"
#include "StridedMatrix.h"

namespace Physica {
    template<class Derived>
    class device_obj<StridedMatrix<Derived>> : public device_obj<LValueMatrix<Derived>> {
        using host_obj = StridedMatrix<Derived>;
        using Base = device_obj<LValueMatrix<Derived>>;
        using This = device_obj<host_obj>;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getRowStride() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getColStride() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getMajorStride() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getMinorStride() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowStrideAtCompile() noexcept { return host_obj::getRowStrideAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getColStrideAtCompile() noexcept { return host_obj::getColStrideAtCompile(); }
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "StridedMatrixImpl/StridedMatrixImpl.cuh"
#include "StridedMatrixImpl/MainDiag.cuh"
#include "StridedMatrixImpl/OffsetDiag.cuh"
