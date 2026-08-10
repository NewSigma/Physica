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

#include "Physica/Macro.h"

struct CUstream_st;
using cudaStream_t = CUstream_st*;

namespace Physica {
    class PHYSICA_API CUDAStream {
        using This = CUDAStream;

        cudaStream_t stream = nullptr;
    public:
        CUDAStream();
        CUDAStream(std::nullptr_t);
        CUDAStream(const This&) = delete;
        CUDAStream(This&& obj) noexcept;
        ~CUDAStream();
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }

        [[nodiscard]] __host__ __device__ operator cudaStream_t() const noexcept;
        /* Operations */
        void wait() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] bool query() const noexcept;
    };
}
