/*
 * Copyright 2024-2026 Weibo He.
 *
 * This file is part of Physica.

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
}
