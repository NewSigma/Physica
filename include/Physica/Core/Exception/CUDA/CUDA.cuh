/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Macro.h"

#ifdef PHYSICA_CUDA
namespace Physica {
    PHYSICA_API [[noreturn]] void throw_cuda_error(int code);

    __host__ __device__ inline void check(cudaError_t err) {
        if (err != cudaSuccess) [[unlikely]] {
        #ifdef __CUDA_ARCH__
            printf("CUDA error with code %d\n", int(err));
        #else
            throw_cuda_error(int(err));
        #endif
        }
    }
}
#endif
