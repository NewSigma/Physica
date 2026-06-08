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

#include <cooperative_groups.h>
#include "Physica/Core/Scalar/Scalar.h"
#include "Physica/Core/Utils/Empty.h"
#include "Physica/Macro.h"

namespace Physica {
    /**
     * \class ThreadBlock is a logical view of cuda thread block. Based on cooperative groups
     */
    template<int NumThread>
    class ThreadBlock final {
        using This = ThreadBlock<NumThread>;
        using Impl = std::conditional_t<NumThread == 1, Empty, cooperative_groups::thread_block>;

        Impl block;
    public:
        __device__ ThreadBlock() requires(NumThread == 1) = default;
        __device__ ThreadBlock() requires(NumThread != 1);
        ThreadBlock(const This&) = default;
        ThreadBlock(This&&) noexcept = default;
        ~ThreadBlock() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        __device__ void sync() const noexcept;
        __device__ auto max(Scalar auto x) const;
        __device__ auto min(Scalar auto x) const;
        __device__ auto sum(Scalar auto x) const;
        /* Static members */
        [[nodiscard]] __device__ constexpr static int tid() noexcept;
        [[nodiscard]] __device__ constexpr static int getNumThread() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getNumThreadAtCompile() noexcept;
    };

    template<int NumThread>
    __device__ inline ThreadBlock<NumThread>::ThreadBlock() requires(NumThread != 1) : block(cooperative_groups::this_thread_block()) {}

    template<int NumThread>
    __device__ void ThreadBlock<NumThread>::sync() const noexcept {
        if constexpr (NumThread != 1)
            block.sync();
    }

    template<int NumThread>
    __device__ auto ThreadBlock<NumThread>::max(Scalar auto x) const {
        using T = decltype(x);
        if constexpr (NumThread == 1)
            return x;
        else {
            static_assert(NumThread != Dynamic, "NoImpl");
            __shared__ std::array<T, NumThread> buffer;
            buffer[tid()] = x;
            const int numThread = getNumThread();
            for (int i = (numThread + 1) / 2; i > 0; i /= 2) {
                if ((tid() < i) && (tid() + i < numThread))
                    buffer[tid()] = std::max(buffer[tid()], buffer[tid() + i]);
                sync();
            }
            return buffer[0];
        }
    }

    template<int NumThread>
    __device__ auto ThreadBlock<NumThread>::min(Scalar auto x) const {
        using T = decltype(x);
        if constexpr (NumThread == 1)
            return x;
        else {
            static_assert(NumThread != Dynamic, "NoImpl");
            __shared__ std::array<T, NumThread> buffer;
            buffer[tid()] = x;
            const int numThread = getNumThread();
            for (int i = (numThread + 1) / 2; i > 0; i /= 2) {
                if ((tid() < i) && (tid() + i < numThread))
                    buffer[tid()] = std::min(buffer[tid()], buffer[tid() + i]);
                sync();
            }
            return buffer[0];
        }
    }

    template<int NumThread>
    __device__ auto ThreadBlock<NumThread>::sum(Scalar auto x) const {
        using T = decltype(x);
        if constexpr (NumThread == 1)
            return x;
        else {
            static_assert(NumThread != Dynamic, "NoImpl");
            __shared__ std::array<T, NumThread> buffer;
            buffer[tid()] = x;
            const int numThread = getNumThread();
            for (int i = (numThread + 1) / 2; i > 0; i /= 2) {
                if ((tid() < i) && (tid() + i < numThread))
                    buffer[tid()] += buffer[tid() + i];
                sync();
            }
            return buffer[0];
        }
    }

    template<int NumThread>
    __device__ constexpr int ThreadBlock<NumThread>::tid() noexcept {
        if constexpr (NumThread != 1)
            return static_cast<int>(Impl::thread_rank());
        return 0;
    }

    template<int NumThread>
    __device__ constexpr int ThreadBlock<NumThread>::getNumThread() noexcept {
        if constexpr (NumThread == Dynamic)
            return static_cast<int>(Impl::num_threads());
        return NumThread;
    }

    template<int NumThread>
    __host__ __device__ consteval int ThreadBlock<NumThread>::getNumThreadAtCompile() noexcept {
        return NumThread;
    }
}
