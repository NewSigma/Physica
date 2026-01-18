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

namespace Physica {
    /**
     * \class ThreadBlock dedicates to
     * 1. Changes logical layout of cuda thread block
     * 2. Wrap cuda cooperative groups
     */
    class ThreadBlock {
        using This = ThreadBlock;
        using Base = cooperative_groups::thread_block;

        Base block;
    public:
        __device__ ThreadBlock();
        ThreadBlock(const This&) = delete;
        ThreadBlock(This&&) noexcept = delete;
        ~ThreadBlock() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ static unsigned int rank() noexcept { return Base::thread_rank(); }
        [[nodiscard]] __device__ static unsigned int getLength() noexcept { return Base::num_threads(); }
    };

    __device__ inline ThreadBlock::ThreadBlock() : block(cooperative_groups::this_thread_block()) {}
}
