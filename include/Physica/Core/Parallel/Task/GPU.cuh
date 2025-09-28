/*
 * Copyright 2025 Weibo He.
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

#include "Physica/Core/Parallel/CUDAContext.cuh"

namespace Physica {
    template<>
    class Task<GPU> {
        using This = Task<GPU>;
    public:
        static void wait();
    };

    inline void Task<GPU>::wait() {
        CUDAContext::getInstance().wait();
    }

    template<ExecutePolicy P>
    Task<Concurrent> parallel_for(std::invocable<size_t> auto fn, size_t num_loop) noexcept requires(P == GPU) {
        return parallel_for<Thread>(fn, num_loop);
    }

    template<ExecutePolicy P>
    Task<Concurrent> parallel_for(std::invocable<size_t> auto fn, size_t num_loop, int part) noexcept requires(P == GPU) {
        return parallel_for<Thread>(fn, num_loop, part);
    }
}
