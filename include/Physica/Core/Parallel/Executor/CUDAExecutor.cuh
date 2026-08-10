/*
 * Copyright 2023-2026 Weibo He.
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

#include "Physica/Core/Utils/NoImpl.h"
#include "Physica/Core/Exception/CUDA/CUDA.cuh"
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "Physica/Core/Parallel/KernelConfig.cuh"

namespace Physica {
    namespace Internal {
        template<class Fn, int MaxThreadsPerBlock, int MinBlocksPerMultiprocessor>
        __global__ void __launch_bounds__(MaxThreadsPerBlock, MinBlocksPerMultiprocessor) kernel(Fn fn) {
            fn();
        }

        template<class Fn, class RtnTy, int MaxThreadsPerBlock, int MinBlocksPerMultiprocessor>
        __global__ void __launch_bounds__(MaxThreadsPerBlock, MinBlocksPerMultiprocessor) kernel(Fn fn, RtnTy* pRet) {
            *pRet = fn();
        }
    }
    /**
     * Single thread with cuda support
     */
    class CUDAExecutor final {
    public:
        template<int MaxThreadsPerBlock = Dynamic, int MinBlocksPerMultiprocessor = Dynamic>
        [[nodiscard]] __host__ __device__ static auto launch(cuda::std::invocable<> auto fn, KernelConfig config, size_t sharedMem = 0);
        [[gnu::nodebug]] static inline void wait();
    private:
        template<int MaxThreadsPerBlock, int MinBlocksPerMultiprocessor>
        __host__ __device__ static void launch_void(cuda::std::invocable<> auto fn, KernelConfig config, size_t sharedMem);

        template<class RtnTy, int MaxThreadsPerBlock, int MinBlocksPerMultiprocessor>
        __host__ __device__ static RtnTy launch_value(cuda::std::invocable<> auto fn, KernelConfig config, size_t sharedMem);
    };

    template<int MaxThreadsPerBlock, int MinBlocksPerMultiprocessor>
    __host__ __device__ auto CUDAExecutor::launch(cuda::std::invocable<> auto fn, KernelConfig config, size_t sharedMem) {
        using RtnTy = cuda::std::invoke_result_t<decltype(fn)>;
        if constexpr (std::same_as<RtnTy, void>)
            return launch_void<MaxThreadsPerBlock, MinBlocksPerMultiprocessor>(std::forward<decltype(fn)>(fn), config, sharedMem);
        else {
            static_assert(std::is_trivially_copyable_v<RtnTy>, "[Error]: Cannot handle non-trivial type");
            return launch_value<RtnTy, MaxThreadsPerBlock, MinBlocksPerMultiprocessor>(std::forward<decltype(fn)>(fn), config, sharedMem);
        }
    }

    template<int MaxThreadsPerBlock, int MinBlocksPerMultiprocessor>
    __host__ __device__ void CUDAExecutor::launch_void(cuda::std::invocable<> auto fn, KernelConfig config, size_t sharedMem) {
        cudaStream_t stream = nullptr;
        // Both host and device must instantiate the global function
        std::ignore = Internal::kernel<std::remove_reference_t<decltype(fn)>, MaxThreadsPerBlock, MinBlocksPerMultiprocessor>;
        if constexpr (IsHost()) {
            stream = cudaStream_t(CUDAContext::getInstance().stream());
            // FIXME: We have to wrap it with lambda because clang rejects it
            [=]() noexcept {
                Internal::kernel<std::remove_reference_t<decltype(fn)>, MaxThreadsPerBlock, MinBlocksPerMultiprocessor><<<config.blocks, config.threads, sharedMem, stream>>>(fn);
            }();
        }
        else
            noImpl("No dynamic parallel support");
        check(cudaGetLastError());
    }
    /**
     * Launch a kernel that yields a trivially copyable type, which is to be returned to the host.
     */
    template<class RtnTy, int MaxThreadsPerBlock, int MinBlocksPerMultiprocessor>
    __host__ __device__ RtnTy CUDAExecutor::launch_value(cuda::std::invocable<> auto fn, KernelConfig config, size_t sharedMem) {
        cudaStream_t stream = nullptr;
        std::ignore = Internal::kernel<std::remove_reference_t<decltype(fn)>, RtnTy, MaxThreadsPerBlock, MinBlocksPerMultiprocessor>;
        if constexpr (IsHost()) {
            stream = cudaStream_t(CUDAContext::getInstance().stream());
            device_obj<Array<RtnTy>> buffer(1);
            [&]() noexcept {
                Internal::kernel<std::remove_reference_t<decltype(fn)>, RtnTy, MaxThreadsPerBlock, MinBlocksPerMultiprocessor><<<config.blocks, config.threads, sharedMem, stream>>>(fn, buffer.data());
            }();
            check(cudaGetLastError());
            return buffer.toHost()[0];
        }
        else
            noImpl("No dynamic parallel support");
    }

    inline void CUDAExecutor::wait() {
        if constexpr (IsHost())
            CUDAContext::getInstance().wait();
    }

    __host__ __device__ inline bool isZeroThread() noexcept {
    #ifdef __CUDA_ARCH__
        return !threadIdx.x && !threadIdx.y && !threadIdx.z && !blockIdx.x && !blockIdx.x && !blockIdx.x;
    #else
        return false;
    #endif
    }
}
