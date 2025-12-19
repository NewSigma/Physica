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
#include "CUDAStream.cuh"
#include "Physica/Core/Exception/CUDA/cuBLAS.cuh"
#include "Physica/Core/Exception/CUDA/cuDSS.cuh"
#include "Physica/Core/Exception/CUDA/cuSolver.cuh"
#include "Physica/Core/Scalar/Scalar.h"

struct cublasContext;

namespace Physica {
    /**
     * \class CUDAContext provides per-thread device resource management
     */
    class PHYSICA_API CUDAContext {
        using This = CUDAContext;

        struct Page {
            int device = 0;
            CUDAStream stream;
            cublasContext* cublas = nullptr;
            cusolverDnContext* cuSolverDn = nullptr;
            cusolverDnParams* cuSolverDnParams = nullptr;
            cudssHandle_t cudss = nullptr;

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
        [[nodiscard]] operator cusolverDnContext*() const noexcept { return pages.top().cuSolverDn; }
        [[nodiscard]] operator cusolverDnParams*() const noexcept { return pages.top().cuSolverDnParams; }
        [[nodiscard]] operator cudssHandle_t() const noexcept { return pages.top().cudss; }
        /* Operations */
        [[nodiscard]] PageGuard push() { return push(device()); }
        [[nodiscard]] PageGuard push(int device);
        [[nodiscard]] cudaError_t query() const noexcept { return getStream().query(); }
        void wait() const { getStream().wait(); }
        /* Getters */
        [[nodiscard]] int device() const noexcept { return pages.top().device; }
        [[nodiscard]] const CUDAStream& getStream() const noexcept { return pages.top().stream; }
        /* Setters */
        void setPointerMode(bool isDeviceSide) const noexcept;
        /* Static members */
        [[nodiscard]] static This& getInstance() noexcept;
        template<Scalar T>
        [[nodiscard]] constexpr static cudaDataType getDataType() noexcept;
    private:
        CUDAContext();
        /* Operations */
        static void setDevice(int device);
    };

    inline void CUDAContext::setPointerMode(bool isDeviceSide) const noexcept {
        cublasSetPointerMode(*this, isDeviceSide ? CUBLAS_POINTER_MODE_DEVICE : CUBLAS_POINTER_MODE_HOST);
    }

    template<Scalar T>
    constexpr cudaDataType CUDAContext::getDataType() noexcept {
        if constexpr (T::isComplex) {
            if constexpr (T::Prec == Float32)
                return cudaDataType::CUDA_C_32F;
            else
                return cudaDataType::CUDA_C_64F;
        }
        else {
            if constexpr (T::Prec == Float32)
                return cudaDataType::CUDA_R_32F;
            else
                return cudaDataType::CUDA_R_64F;
        }
    }

    __host__ __device__ inline bool isZeroThread() {
    #ifdef __CUDA_ARCH__
        return !threadIdx.x && !threadIdx.y && !threadIdx.z && !blockIdx.x && !blockIdx.x && !blockIdx.x;
    #else
        return false;
    #endif
    }
}
