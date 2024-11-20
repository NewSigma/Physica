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

#include "Physica/Core/Exception/CUDA/cuBLAS.cuh"
#include "CUDAStream.cuh"

struct cublasContext;

namespace Physica::Core {
    /**
     * \class CUDAContext provides per-thread device resource management
     */
    class PHYSICA_API CUDAContext : private CUDAStream {
        using This = CUDAContext;
        using Base = CUDAStream;

        cublasContext* cublas;
    public:
        CUDAContext(const This&) = delete;
        CUDAContext(This&&) noexcept = delete;
        ~CUDAContext();
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] operator cudaStream_t() const noexcept { return getStream(); }
        [[nodiscard]] operator cublasContext*();
        /* Operations */
        using Base::query;
        using Base::wait;
        /* Getters */
        [[nodiscard]] const CUDAStream& getStream() const noexcept { return *this; }
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
}
