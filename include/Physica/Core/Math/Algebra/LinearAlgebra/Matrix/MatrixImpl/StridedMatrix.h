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

    template<class Derived>
    constexpr size_t StridedMatrix<Derived>::getRowStride() const noexcept {
        if constexpr (Derived::getRowStrideAtCompile() != Dynamic)
            return Derived::getRowStrideAtCompile();
        else
            return Base::getDerived().getRowStride();
    }

    template<class Derived>
    constexpr size_t StridedMatrix<Derived>::getColStride() const noexcept {
        if constexpr (Derived::getColStrideAtCompile() != Dynamic)
            return Derived::getColStrideAtCompile();
        else
            return Base::getDerived().getColStride();
    }

    template<class Derived>
    auto StridedMatrix<Derived>::data_handle() noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    auto StridedMatrix<Derived>::data_handle() const noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    auto StridedMatrix<Derived>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        return self.data_handle() + row * self.getRowStride() + col * self.getColStride();
    }

    template<class Derived>
    StridedMatrix<Derived>::StridedMatrix() {
        static_assert(Derived::isRowMatrix() || Derived::isColMatrix(), "[Error]: StridedMatrix requires a known major");
    }
}
