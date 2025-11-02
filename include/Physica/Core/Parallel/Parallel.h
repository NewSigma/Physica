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

#include "Physica/Macro.h"

namespace Physica {
    enum ExecutePolicy : char {
        Sequential,
        Concurrent,
        Thread,
        GPU
    };
    /**
     * \class Task maintains lifetime of async tasks
     */
    template<ExecutePolicy P>
    class [[nodiscard]] Task;

    template<ExecutePolicy P>
    [[nodiscard]] Task<P> schedule(std::invocable<> auto fn) noexcept {
        static_assert(P != Concurrent, "[Error]: Not support");
        fn();
        co_return;
    }

    template<ExecutePolicy P>
    void parallel_for(auto, size_t) noexcept {
        static_assert(P == Sequential, "[Error]: Not implemented");
    }

    template<ExecutePolicy P>
    void parallel_for(auto, size_t, int) noexcept {
        static_assert(P == Sequential, "[Error]: Not implemented");
    }
}

#include "Task/Sequential.h"
#include "Task/Concurrent.h"
#include "Task/Thread.h"
#ifdef PHYSICA_CUDA
    #include "Task/GPU.cuh"
#endif
