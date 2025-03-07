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

#include <stack>
#include "Physica/Core/Exception/CUDA/cuBLAS.cuh"
#include "CUDAStream.cuh"

struct cublasContext;

namespace Physica {
    /**
     * \class CUDAContext provides per-thread device resource management
     */
    class PHYSICA_API CUDAContext {
        using This = CUDAContext;

        struct Page {
            int device;
            CUDAStream stream;
            cublasContext* cublas;

            Page(int device_, CUDAStream stream_);
            Page(const Page&) = delete;
            Page(Page&&) noexcept = default;
            ~Page();
            /* Operators */
            Page& operator=(const Page&) = delete;
            Page& operator=(Page&&) noexcept = delete;
        };

        struct PageGuard {
            ~PageGuard() { CUDAContext::getInstance().pages.pop(); }
        };

        std::stack<Page> pages;
    public:
        CUDAContext(const This&) = delete;
        CUDAContext(This&&) noexcept = delete;
        ~CUDAContext();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] operator cudaStream_t() const noexcept { return getStream(); }
        [[nodiscard]] operator cublasContext*() const noexcept { return pages.top().cublas; }
        /* Operations */
        PageGuard push() { return push(device()); }
        PageGuard push(int device);
        [[nodiscard]] cudaError_t query() { return getStream().query(); }
        void wait() { getStream().wait(); }
        /* Getters */
        [[nodiscard]] int device() const noexcept { return pages.top().device; }
        [[nodiscard]] const CUDAStream& getStream() const noexcept { return pages.top().stream; }
        /* Setters */
        void setPointerMode(bool isDeviceSide) noexcept;
        /* Static members */
        [[nodiscard]] static This& getInstance() noexcept;
    private:
        CUDAContext();
    };

    inline void CUDAContext::setPointerMode(bool isDeviceSide) noexcept {
        cublasSetPointerMode(*this, isDeviceSide ? CUBLAS_POINTER_MODE_DEVICE : CUBLAS_POINTER_MODE_HOST);
    }

    __host__ __device__ inline bool isZeroThread() {
    #ifdef __CUDA_ARCH__
        return !threadIdx.x && !threadIdx.y && !threadIdx.z && !blockIdx.x && !blockIdx.x && !blockIdx.x;
    #else
        return false;
    #endif
    }
}
