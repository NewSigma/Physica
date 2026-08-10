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

#include <stack>
#include "Physica/Core/Parallel/CUDAStream.cuh"
#include "Physica/Core/Scalar/Scalar.h"

struct cublasContext;
struct cusolverDnContext;
struct cusolverDnParams;
struct cudssContext;

namespace Physica {
    /**
     * \class CUDAContext provides per-thread device resource management
     */
    class PHYSICA_API CUDAContext {
        using This = CUDAContext;
        struct Page {
            int device = 0;
        #ifdef PHYSICA_CUDA
            CUDAStream stream;
            cublasContext* cublas = nullptr;
            cusolverDnContext* cuSolverDn = nullptr;
            cusolverDnParams* cuSolverDnParams = nullptr;
            cudssContext* cudss = nullptr;
        #endif

            Page(int device_, CUDAStream stream_);
            Page(const Page&) = delete;
            Page(Page&&) noexcept = default;
            ~Page();
            /* Operators */
            Page& operator=(const Page&) = delete;
            Page& operator=(Page&&) noexcept = delete;
        };

        std::stack<Page> pages;
    public:
        struct PageGuard {
            ~PageGuard() { CUDAContext::getInstance().pages.pop(); }
        };
    public:
        CUDAContext(const This&) = delete;
        CUDAContext(This&&) noexcept = delete;
        ~CUDAContext();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] operator cudaStream_t() const noexcept;
        [[nodiscard]] operator cublasContext*() const noexcept;
        [[nodiscard]] operator cusolverDnContext*() const noexcept;
        [[nodiscard]] operator cusolverDnParams*() const noexcept;
        [[nodiscard]] operator cudssContext*() const noexcept;
        /* Operations */
        [[nodiscard]] PageGuard push();
        [[nodiscard]] PageGuard push(int device);
        [[gnu::nodebug]] void wait() const;
        /* Getters */
        [[nodiscard]] int device() const noexcept;
        [[nodiscard]] const CUDAStream& stream() const noexcept;
        /* Setters */
        void setPointerMode(bool isDeviceSide) const noexcept;
        /* Static members */
        [[nodiscard]] static This& getInstance() noexcept;
        template<Scalar T>
        [[nodiscard]] constexpr static auto getDataType() noexcept;
    private:
        CUDAContext();
        /* Operations */
        static void setDevice(int device);
    };

    template<Scalar T>
    constexpr auto CUDAContext::getDataType() noexcept {
    #ifdef PHYSICA_CUDA
        if constexpr (T::isComplex()) {
            if constexpr (T::Prec == Float32)
                return CUDA_C_32F;
            else
                return CUDA_C_64F;
        }
        else {
            if constexpr (T::Prec == Float32)
                return CUDA_R_32F;
            else
                return CUDA_R_64F;
        }
    #else
        static_assert(false, "[Error]: No impl");
    #endif
    }
}
