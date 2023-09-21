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
#include <cassert>
#include <thread>
#include "Physica/Core/Parallel/CudaStream.cuh"
#include "Physica/Utils/CUDA/DebugUtil.cuh"

namespace Physica::Core {
    CudaStream::CudaStream() {
        cudaStreamCreate(&stream);
    }

    CudaStream::CudaStream(cudaStream_t stream_) : stream(stream_) {}

    CudaStream::CudaStream(CudaStream&& obj) noexcept : stream(obj.stream) {
        obj.stream = nullptr;
    }

    CudaStream::~CudaStream() {
        cudaStreamDestroy(stream);
        stream = nullptr;
    }

    CudaStream& CudaStream::operator=(CudaStream obj) noexcept {
        swap(obj);
        return *this;
    }

    cudaError_t CudaStream::query() const {
        return cudaStreamQuery(stream);
    }

    void CudaStream::wait() const {
        while (query() != cudaSuccess)
            std::this_thread::yield();
        cudaCheck(cudaGetLastError());
    }

    void CudaStream::swap(CudaStream& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(stream, obj.stream);
    }

    CudaStream CudaStream::makeStream() {
        cudaStream_t stream = nullptr;
        cudaStreamCreate(&stream);
        return CudaStream(stream);
    }
}
