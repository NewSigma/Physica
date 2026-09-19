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

#include "RValueTensor.cuh"
#include "LValueTensor.h"

namespace Physica {
    template<class Derived>
    class device_obj<LValueTensor<Derived>> : public device_obj<RValueTensor<Derived>> {
        using This = device_obj<LValueTensor<Derived>>;
        using Base = device_obj<RValueTensor<Derived>>;
    public:
        using typename Base::IndexType;
    protected:
        using typename Base::T;
        using typename Base::Trv;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;

        __host__ __device__ auto operator=(Scalar auto x) -> device_obj<Derived>&;
        __host__ __device__ void operator+=(Scalar auto x);
        __host__ __device__ void operator-=(Scalar auto x);
        __host__ __device__ void operator*=(Scalar auto x);
        __host__ __device__ void operator/=(Scalar auto x);

        __host__ __device__ auto operator=(const Tensor auto& x) -> device_obj<Derived>&;
        __host__ __device__ void operator+=(const Tensor auto& x);
        __host__ __device__ void operator-=(const Tensor auto& x);

        [[nodiscard]] __device__ decltype(auto) operator[](this auto&&, const IndexType& index);
        [[nodiscard]] __device__ decltype(auto) operator[](this auto&&, std::integral auto... dims);
        /* Operations */
        [[nodiscard]] __device__ decltype(auto) calc(const IndexType& index) const { return operator[](index); }
        [[nodiscard]] __device__ decltype(auto) calc(std::integral auto... dims) const { return operator[](dims...); }

        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        /* Getters */
        [[nodiscard]] __device__ auto data_ptr(this auto&&, const IndexType& index) noexcept;
        [[nodiscard]] __device__ auto data_ptr(this auto&&, std::integral auto... dims) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isLValueTensor() noexcept { return true; }
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "LValueTensorImpl/LValueTensorImpl.cuh"
