/*
 * Copyright 2023-2024 Weibo He.
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

#include <Physica/Core/Parallel/CUDAContext.cuh>
#include <Physica/Core/Exception/CUDA/CUDA.cuh>
#include "SequentialExecutor.h"

namespace Physica::Core {
    /**
     * Single thread with cuda support
     */
    class CUDAExecutor : public SequentialExecutor {
    public:
        static void wait() { check(cudaDeviceSynchronize()); }
    };
}

namespace Physica {
    template<>
    class Traits<Core::CUDAExecutor> {
    public:
        constexpr static bool UseCPU = false;
        constexpr static bool UseCUDA = true;
    };
}
