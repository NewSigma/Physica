/*
 * Copyright 2024-2025 Weibo He.
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
    __host__ __device__ inline Real<Float16> abs(const Real<Float16>& x) noexcept {
        return Real<Float16>(::__habs(x.toMachine()));
    }

    template<>
    __host__ __device__ inline Real<Float16> sqrt(const Real<Float16>& x) noexcept {
        if constexpr (IsHost())
            return float16(sqrt(float32(x)));
        else {
        #ifdef __CUDA_ARCH__
            return Real<Float16>(::hsqrt(x.toMachine()));
        #endif
        }
    }

    template<>
    __host__ __device__ inline Real<Float16> ln(const Real<Float16>& x) noexcept {
        if constexpr (IsHost())
            return float16(ln(float32(x)));
        else {
        #ifdef __CUDA_ARCH__
            return Real<Float16>(::hlog(x.toMachine()));
        #endif
        }
    }

    template<>
    __host__ __device__ inline Real<Float16> exp(const Real<Float16>& x) noexcept {
        if constexpr (IsHost())
            return float16(exp(float32(x)));
        else {
        #ifdef __CUDA_ARCH__
            return Real<Float16>(::hexp(x.toMachine()));
        #endif
        }
    }
}
