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
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "Physica/Core/Exception/CUDA/cuBLAS.cuh"

namespace Physica::Core {
    CUDAContext::CUDAContext() : CUDAStream(), cublas(nullptr) {}

    CUDAContext::~CUDAContext() {
        cublasDestroy(cublas);
    }

    CUDAContext::operator cublasContext*() {
        if (cublas == nullptr) {
            check(cublasCreate(&cublas));
            check(cublasSetStream(cublas, getStream()));
        }
        return cublas;
    }

    CUDAContext& CUDAContext::getInstance() noexcept {
        thread_local static auto* instance = new CUDAContext();
        return *instance;
    }
}
