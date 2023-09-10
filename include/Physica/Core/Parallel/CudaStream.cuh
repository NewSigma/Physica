/*
 * Copyright 2023 WeiBo He.
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

namespace Physica::Core {
    class StreamPool;

    class CudaStream {
        cudaStream_t stream;
    public:
        CudaStream();
        CudaStream(const CudaStream&) = delete;
        CudaStream(CudaStream&& obj) noexcept;
        ~CudaStream();
        /* Operators */
        CudaStream& operator=(CudaStream obj) noexcept;
        [[nodiscard]] operator cudaStream_t() const noexcept { return stream; }
        /* Operations */
        [[nodiscard]] cudaError_t query() const;
        void wait() const;
        void swap(CudaStream& obj) noexcept;
        /* Static members */
        [[nodiscard]] static CudaStream makeStream();
    private:
        CudaStream(cudaStream_t stream_);

        friend class StreamPool;
    };
}
