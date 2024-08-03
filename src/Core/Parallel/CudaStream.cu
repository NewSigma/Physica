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
#include <cassert>
#include <thread>
#include <Physica/Core/Parallel/CudaStream.cuh>
#include <Physica/Core/Exception/CudaException.cuh>

namespace Physica::Core {
    CudaStream::CudaStream() {
        cudaStreamCreate(&stream);
    }

    CudaStream::CudaStream(CudaStream&& obj) noexcept : stream(obj.stream) {
        obj.stream = nullptr;
    }

    CudaStream::~CudaStream() {
        cudaStreamDestroy(stream);
        stream = nullptr;
    }

    cudaError_t CudaStream::query() const {
        return cudaStreamQuery(stream);
    }

    void CudaStream::wait() const {
        cudaCheck(cudaStreamSynchronize(stream));
    }

    void CudaStream::swap(CudaStream& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(stream, obj.stream);
    }
}
