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

#include "LValueMatrix.h"

namespace Physica {
    template<class Derived>
    class StridedMatrix : public LValueMatrix<Derived> {
        using Base = LValueMatrix<Derived>;
        using This = StridedMatrix<Derived>;
    public:
        ~StridedMatrix() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;
        /* Getters */
        [[nodiscard]] constexpr size_t getRowStride() const noexcept;
        [[nodiscard]] constexpr size_t getColStride() const noexcept;
        [[nodiscard]] constexpr size_t getMajorStride() const noexcept;
        [[nodiscard]] constexpr size_t getMinorStride() const noexcept;
        [[nodiscard]] auto data_handle() noexcept;
        [[nodiscard]] auto data_handle() const noexcept;
        [[nodiscard]] auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStrided() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static size_t getRowStrideAtCompile() noexcept { return Dynamic; }
        [[nodiscard]] __host__ __device__ consteval static size_t getColStrideAtCompile() noexcept { return Dynamic; }
    protected:
        StridedMatrix();
        StridedMatrix(const This&) = default;
        StridedMatrix(This&&) noexcept = default;
    };
}

#include "StridedMatrixImpl/StridedMatrixImpl.h"
#include "StridedMatrixImpl/MainDiag.h"
#include "StridedMatrixImpl/OffsetDiag.h"
