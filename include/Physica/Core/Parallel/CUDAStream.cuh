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

#include <cuda_runtime.h>
#include "Physica/Macro.h"

namespace Physica::Core {
    class PHYSICA_API CUDAStream {
        cudaStream_t stream;
    public:
        CUDAStream();
        CUDAStream(const CUDAStream&) = delete;
        CUDAStream(CUDAStream&& obj) noexcept;
        ~CUDAStream();
        /* Operators */
        CUDAStream& operator=(CUDAStream obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] operator cudaStream_t() const noexcept { return stream; }
        /* Operations */
        [[nodiscard]] cudaError_t query() const;
        void wait() const;
        void swap(CUDAStream& __restrict obj) noexcept;
    };
}
