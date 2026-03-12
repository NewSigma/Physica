/*
 * Copyright 2024-2026 Weibo He.
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
#include <cstring>
#include "Physica/Macro.h"

namespace Physica {
    namespace Internal {
        struct Undef final {
            template<class T>
            __host__ __device__ constexpr operator T() const noexcept {
                std::array<char, sizeof(T)> bytes;
                return std::bit_cast<T>(bytes);
            }
        };
    }

    [[gnu::const, gnu::nodebug]] __host__ __device__ constexpr auto undef() noexcept {
        return Internal::Undef{};
    }
    /**
     * Improvements:
     * 1. Emit assertion to help debugging
     * 2. __host__ __device__ make it available in kernel functions
     */
    [[noreturn, gnu::always_inline]] __host__ __device__ inline void unreachable([[maybe_unused]] const char* msg = nullptr) noexcept {
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

    [[gnu::always_inline, gnu::nodebug]] __host__ __device__ inline void assume(bool predicate) noexcept {
        assert(predicate && "[Error]: Bad assumption");
        [[assume(predicate)]];
    }

    template<class T>
    void memswap(T* a, T* b) noexcept {
        // FIXME: static_assert trivially_relocatable once we dump to CXX26
        assert(a != b && "[Error]: Self swap is likely a bug");
        alignas(T) std::array<std::byte, sizeof(T)> buffer;
        memcpy(buffer.data(), (void*)a, buffer.size());
        memcpy((void*)a, (void*)b, buffer.size());
        memcpy((void*)b, buffer.data(), buffer.size());
    }
}
