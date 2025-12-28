/*
 * Copyright 2024-2025 Weibo He.
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

#include <cassert>
#include "Physica/Macro.h"

namespace Physica {
    /**
     * Improvements:
     * 1. Emit assertion to help debugging
     * 2. __host__ __device__ make it available in kernel functions
     */
    [[noreturn, gnu::always_inline]] __host__ __device__ inline void unreachable() noexcept {
        assert(false && "[Error]: Trigger unreachable");
    #ifdef __CUDA_ARCH__
        __builtin_unreachable();
    #else
        if constexpr (IsHost()) {
        #if defined(_MSC_VER) && !defined(__clang__)
            __assume(false);
        #else
            __builtin_unreachable();
        #endif
        }
    #endif
    }
}
