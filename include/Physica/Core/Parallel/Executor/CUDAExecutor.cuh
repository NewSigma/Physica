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
#include "SeqExecutor.h"

namespace Physica {
    namespace Internal {
        template<class Functor>
        __global__ void kernel(Functor func) {
            return func();
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
    class CUDAExecutor : public SeqExecutor {
        struct KernelFuture {
            cudaStream_t stream;

            void wait() { check(cudaStreamSynchronize(stream)); }
        };
    public:
        template<class Functor>
        __host__ __device__ static KernelFuture launch(Functor func, KernelConfig config, size_t sharedMem = 0);
        static inline void wait();
    };

    template<class Functor>
    __host__ __device__ auto CUDAExecutor::launch(Functor func, KernelConfig config, size_t sharedMem) -> KernelFuture {
        cudaStream_t stream = nullptr;
        if constexpr (IsHost())
            stream = cudaStream_t(CUDAContext::getInstance().getStream());
        else
            noImpl("No dynamic parallel support");
        Internal::kernel<<<config.blocks, config.threads, sharedMem, stream>>>(func);
        check(cudaGetLastError());
        return KernelFuture(stream);
    }

    inline void CUDAExecutor::wait() {
        if constexpr (IsHost())
            CUDAContext::getInstance().wait();
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

namespace std {
    template<>
    struct formatter<Physica::KernelConfig, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        auto format(const Physica::KernelConfig& obj, std::format_context& ctx) const {
            return std::format_to(ctx.out(), "Blocks: ({}, {}, {})\nThreads: ({}, {}, {})",
                    obj.blocks.x, obj.blocks.y, obj.blocks.z, obj.threads.x, obj.threads.y, obj.threads.z);
        }
    };
}
