/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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

__host__ __device__ CUDAStream::operator cudaStream_t() const noexcept {
    return reinterpret_cast<cudaStream_t>(stream);
}

bool CUDAStream::query() const noexcept {
    assert(stream != nullptr);
    return cudaStreamQuery(stream) == cudaSuccess;
}

void CUDAStream::wait() const {
    check(cudaStreamSynchronize(stream));
}

void CUDAStream::swap(CUDAStream& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    std::swap(stream, obj.stream);
}
