/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Half2.h"

namespace Physica {
    [[nodiscard]] __host__ __device__ inline SIMD<Real<Float16>, 2> fma(
            const SIMD<Real<Float16>, 2> a, const SIMD<Real<Float16>, 2> b, const SIMD<Real<Float16>, 2> c) noexcept {
    #ifdef __CUDA_ARCH__
        return SIMD<Real<Float16>, 2>(__hfma2(a.toMachine(), b.toMachine(), c.toMachine()));
    #else
        return a * b + c;
    #endif
    }
}
