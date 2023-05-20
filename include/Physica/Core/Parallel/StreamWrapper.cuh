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

    class StreamWrapper {
        cudaStream_t stream;
    public:
        StreamWrapper();
        StreamWrapper(const StreamWrapper&) = delete;
        StreamWrapper(StreamWrapper&& obj) noexcept;
        ~StreamWrapper();
        /* Operators */
        StreamWrapper& operator=(StreamWrapper obj) noexcept;
        /* Operations */
        void swap(StreamWrapper& obj) noexcept;
        /* Getters */
        [[nodiscard]] cudaStream_t getStream() const noexcept { return stream; }
    private:
        StreamWrapper(cudaStream_t stream_);

        friend class StreamPool;
    };
}
