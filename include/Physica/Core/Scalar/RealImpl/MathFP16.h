/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Float16.h"

namespace Physica {
    template<>
    __host__ __device__ inline float16 abs(const float16& x) noexcept {
        return float16(::__habs(x.toMachine()));
    }

    template<>
    __host__ __device__ inline float16 sqrt(const float16& x) noexcept {
        if constexpr (IsHost())
            return float16(sqrt(float32(x)));
        else {
        #ifdef __CUDA_ARCH__
            return float16(::hsqrt(x.toMachine()));
        #else
            unreachable();
        #endif
        }
    }

    template<>
    __host__ __device__ inline float16 ln(const float16& x) noexcept {
        if constexpr (IsHost())
            return float16(ln(float32(x)));
        else {
        #ifdef __CUDA_ARCH__
            return float16(::hlog(x.toMachine()));
        #else
            unreachable();
        #endif
        }
    }

    template<>
    __host__ __device__ inline float16 ln1p(const float16& x) noexcept {
        return float16(ln1p(float32(x)));
    }

    template<>
    __host__ __device__ inline float16 exp(const float16& x) noexcept {
        if constexpr (IsHost())
            return float16(exp(float32(x)));
        else {
        #ifdef __CUDA_ARCH__
            return float16(::hexp(x.toMachine()));
        #else
            unreachable();
        #endif
        }
    }

    template<>
    __host__ __device__ inline float16 sinh(const float16& x) noexcept {
        return float16(sinh(float32(x)));
    }
}
