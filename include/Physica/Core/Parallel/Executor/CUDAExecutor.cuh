/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/Core/Exception/CUDA/CUDA.cuh"
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "SeqExecutor.h"

namespace Physica {
    namespace Internal {
        template<class Functor>
        __global__ void kernel(Functor func) {
            return func();
        }
    }
    /**
     * Single thread with cuda support
     */
    class CUDAExecutor : public SeqExecutor {
    public:
        template<class Functor>
        static CUDAExecutor launch(Functor func, dim3 numBlocks, dim3 numThreads, size_t sharedMem = 0);

        static void wait() {
        #ifdef PHYSICA_CUDA
            check(cudaDeviceSynchronize());
        #endif
        }
    };

    template<class Functor>
    CUDAExecutor CUDAExecutor::launch(Functor func, dim3 numBlocks, dim3 numThreads, size_t sharedMem) {
        Internal::kernel<<<numBlocks, numThreads, sharedMem, CUDAContext::getInstance()>>>(func);
        check(cudaGetLastError());
        return {};
    }
}

namespace Physica {
    template<>
    class Traits<CUDAExecutor> {
    public:
        constexpr static bool UseCPU = false;
        constexpr static bool UseCUDA = true;
    };
}
