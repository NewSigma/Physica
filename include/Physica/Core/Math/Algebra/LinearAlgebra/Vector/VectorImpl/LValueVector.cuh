/*
 * Copyright 2022-2026 Weibo He.
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

#include "LValueVector.h"
#include "RValueVector.cuh"

namespace Physica {
    template<class Derived>
    class device_obj<LValueVector<Derived>> : public device_obj<RValueVector<Derived>> {
        using This = device_obj<LValueVector<Derived>>;
        using Base = device_obj<RValueVector<Derived>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) = delete;

        __host__ __device__ auto operator=(Scalar auto x) -> device_obj<Derived>&;
        __host__ __device__ void operator+=(Scalar auto x);
        __host__ __device__ void operator-=(Scalar auto x);
        __host__ __device__ void operator*=(Scalar auto x);
        __host__ __device__ void operator/=(Scalar auto x);

        __host__ __device__ device_obj<Derived>& operator=(const Vector auto& v);
        __host__ __device__ void operator+=(const Vector auto& v);
        __host__ __device__ void operator-=(const Vector auto& v);

        [[nodiscard]] __device__ decltype(auto) operator[](this auto&&, size_t index);
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t index) const;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto) const;

        __host__ __device__ void reverse(const auto& grad) const noexcept;

        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto head(this auto&&, size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto tail(this auto&&, size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] __host__ __device__ auto segment(this auto&&, size_t from, size_t to) noexcept;
        template<int Major, size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto reshape(this auto&& self, size_t row, size_t col) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto reshape_row(this auto&& self, size_t row, size_t col) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto reshape_col(this auto&& self, size_t row, size_t col) noexcept;
        [[nodiscard]] __host__ __device__ auto reshape_like(this auto&& self, const Matrix auto& mat) noexcept;

        void zero_grad() noexcept;
        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_any(auto& distribution);
        /* Getters */
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t index) noexcept;
        [[nodiscard]] __host__ __device__ auto& front(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto& back(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isLValueVector() noexcept { return true; }
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "LValueVectorImpl/LValueVectorImpl.cuh"
#include "LValueVectorImpl/VectorConvert/RealVector.cuh"
#include "LValueVectorImpl/VectorConvert/ImagVector.cuh"
#include "LValueVectorImpl/LVectorBlock.cuh"
