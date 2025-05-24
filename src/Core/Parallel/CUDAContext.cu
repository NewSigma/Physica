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
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "Physica/Core/Exception/CUDA/CUDA.cuh"
#include "Physica/Core/Exception/CUDA/cuBLAS.cuh"

namespace Physica {
    CUDAContext::Page::Page(int device_, CUDAStream stream_) : device(device_), stream(std::move(stream_)) {
        check(cublasCreate(&cublas));
        check(cublasSetStream(cublas, stream));

        check(cusolverDnCreate(&cuSolverDn));
        check(cusolverDnSetStream(cuSolverDn, stream));

        check(cusolverDnCreateParams(&cuSolverDnParams));
    }

    CUDAContext::Page::~Page() {
        cublasDestroy(cublas);
        cusolverDnDestroy(cuSolverDn);
        cusolverDnDestroyParams(cuSolverDnParams);
    }

    CUDAContext::CUDAContext() {
        setDevice(0);
        pages.emplace(0, CUDAStream(nullptr));
    }

    CUDAContext::~CUDAContext() {
        while (!pages.empty())
            pages.pop();
    }

    auto CUDAContext::push(int device) -> PageGuard {
        setDevice(0);
        pages.emplace(device, CUDAStream());
        return {};
    }

    CUDAContext& CUDAContext::getInstance() noexcept {
        thread_local static auto* instance = new CUDAContext();
        return *instance;
    }

    void CUDAContext::setDevice(int device) {
        check(cudaSetDevice(device));

        uint64_t val = std::numeric_limits<uint64_t>::max();
        cudaMemPool_t memPool;
        check(cudaDeviceGetDefaultMemPool(&memPool, device));
        check(cudaMemPoolSetAttribute(memPool, cudaMemPoolAttrReleaseThreshold, &val));
    }
}
