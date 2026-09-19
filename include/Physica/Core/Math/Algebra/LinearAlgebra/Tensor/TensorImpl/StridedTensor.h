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

#include "LValueTensor.h"

namespace Physica {
    template<class Derived>
    class StridedTensor : public LValueTensor<Derived> {
        using Base = LValueTensor<Derived>;
        using This = StridedTensor<Derived>;
    public:
        using typename Base::IndexType;
    public:
        ~StridedTensor() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;
        /* Getters */
        [[nodiscard]] auto getStrides() const noexcept;
        [[nodiscard]] size_t getStride(size_t dim) const noexcept;
        [[nodiscard]] auto data_handle() noexcept;
        [[nodiscard]] auto data_handle() const noexcept;
        [[nodiscard]] auto data_ptr(this auto&&, const IndexType& index) noexcept;
        [[nodiscard]] auto data_ptr(this auto&&, std::integral auto... dims) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStrided() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static IndexType getStrideAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getStrideAtCompile(size_t dim) noexcept;
    protected:
        StridedTensor() = default;
        StridedTensor(const This&) = default;
        StridedTensor(This&&) noexcept = default;
    };
}

#include "StridedTensorImpl/StridedTensorImpl.h"
