/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "StridedTensor.cuh"
#include "CompactTensor.h"

namespace Physica {
    template<class Derived>
    class device_obj<CompactTensor<Derived>> : public device_obj<StridedTensor<Derived>> {
        using host_obj = CompactTensor<Derived>;
        using This = device_obj<host_obj>;
        using Base = device_obj<StridedTensor<Derived>>;
    public:
        using typename Base::IndexType;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        __host__ __device__ void assign(this const auto& self, Tensor auto&& target);
        /* Getters */
        [[nodiscard]] __host__ __device__ auto data(this auto&& self) noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr auto getStrides() const noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "CompactTensorImpl/CompactTensorImpl.cuh"
#include "CompactTensorImpl/Flatten.cuh"
#include "CompactTensorImpl/TensorSlice.cuh"
#include "CompactTensorImpl/TensorBlock.cuh"
