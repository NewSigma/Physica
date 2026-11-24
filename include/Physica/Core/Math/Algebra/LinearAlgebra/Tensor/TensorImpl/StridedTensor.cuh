/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
    protected:
        using typename Base::T;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;
        /* Operations */
        __host__ __device__ void assign(this const auto& self, Tensor auto&& target);
        __host__ __device__ void zeros() noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto getStrides() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getStride(size_t dim) const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, const IndexType& index) noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, std::same_as<size_t> auto... dims) noexcept;
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
#include "StridedTensorImpl/TensorFiber.cuh"
