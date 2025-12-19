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
#include <cassert>
#include "Physica/Core/Exception/CUDA/CUDA.cuh"
#include "Physica/Core/Parallel/CUDAStream.cuh"

using namespace Physica;

CUDAStream::CUDAStream() {
    cudaStreamCreate(&stream);
}

CUDAStream::CUDAStream(std::nullptr_t) {}

CUDAStream::CUDAStream(CUDAStream&& obj) noexcept : stream(obj.stream) {
    obj.stream = nullptr;
}

CUDAStream::~CUDAStream() {
    if (stream != nullptr) {
        cudaStreamDestroy(stream);
        stream = nullptr;
    }
}

cudaError_t CUDAStream::query() const noexcept {
    return cudaStreamQuery(stream);
}

void CUDAStream::wait() const {
    check(cudaStreamSynchronize(stream));
}

void CUDAStream::swap(CUDAStream& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    std::swap(stream, obj.stream);
}
