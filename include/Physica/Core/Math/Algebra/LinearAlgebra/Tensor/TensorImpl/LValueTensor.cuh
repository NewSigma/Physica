/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
        using typename Base::Tv;
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
        [[nodiscard]] __device__ decltype(auto) operator[](this auto&&, std::same_as<size_t> auto... dims);
        /* Operations */
        [[nodiscard]] __device__ decltype(auto) calc(const IndexType& index) const { return operator[](index); }
        [[nodiscard]] __device__ decltype(auto) calc(std::same_as<size_t> auto... dims) const { return operator[](dims...); }

        __device__ void forND(std::invocable<T&, IndexType> auto fn);
        __device__ void forND(std::invocable<const T&, IndexType> auto fn) const;

        __host__ __device__ void zero_grad() noexcept;
        __host__ __device__ void zeros() noexcept;
        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_any(auto& distribution);
        /* Getters */
        [[nodiscard]] __device__ auto data_ptr(this auto&&, const IndexType& index) noexcept;
        [[nodiscard]] __device__ auto data_ptr(this auto&&, std::same_as<size_t> auto... dims) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isLValueTensor() noexcept { return true; }
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "LValueTensorImpl/LValueTensorImpl.cuh"
#include "LValueTensorImpl/Flatten.cuh"
#include "LValueTensorImpl/TensorFiber.cuh"
#include "LValueTensorImpl/TensorSlice.cuh"
#include "LValueTensorImpl/TensorBlock.cuh"
