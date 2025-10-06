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

#include <format>
#include "Physica/Core/Exception/NoImplException.h"
#include "Physica/Core/Exception/CUDA/CUDA.cuh"
#include "Physica/Core/Parallel/CUDAContext.cuh"

namespace Physica {
    namespace Internal {
        template<int MaxThreadsPerBlock, int MinBlocksPerMultiprocessor>
        __global__ void __launch_bounds__(MaxThreadsPerBlock, MinBlocksPerMultiprocessor) kernel(auto fn) {
            return fn();
        }
    }

    struct KernelConfig {
        dim3 blocks;
        dim3 threads;

        __host__ __device__ KernelConfig(dim3 blocks_, dim3 threads_) : blocks(blocks_), threads(threads_) {}
    };
    /**
     * Single thread with cuda support
     */
    class CUDAExecutor {
        struct KernelFuture {
            cudaStream_t stream;

            void wait() { check(cudaStreamSynchronize(stream)); }
        };
    public:
        template<int MaxThreadsPerBlock = Dynamic, int MinBlocksPerMultiprocessor = Dynamic>
        __host__ __device__ static KernelFuture launch(auto&& fn, KernelConfig config, size_t sharedMem = 0);
        static inline void wait();
    };

    template<int MaxThreadsPerBlock, int MinBlocksPerMultiprocessor>
    __host__ __device__ auto CUDAExecutor::launch(auto&& fn, KernelConfig config, size_t sharedMem) -> KernelFuture {
        cudaStream_t stream = nullptr;
        if constexpr (IsHost())
            stream = cudaStream_t(CUDAContext::getInstance().getStream());
        else
            noImpl("No dynamic parallel support");
        Internal::kernel<MaxThreadsPerBlock, MinBlocksPerMultiprocessor><<<config.blocks, config.threads, sharedMem, stream>>>(fn);
        check(cudaGetLastError());
        return KernelFuture(stream);
    }

    inline void CUDAExecutor::wait() {
        if constexpr (IsHost())
            CUDAContext::getInstance().wait();
    }
}

namespace std {
    template<>
    struct PHYSICA_API formatter<Physica::KernelConfig, char> {
        constexpr static auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
        static auto format(const Physica::KernelConfig& obj, std::format_context& ctx) -> std::format_context::iterator;
    };
}
