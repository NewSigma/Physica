/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../SIMD.cuh"

namespace Physica {
    template<FloatPrec Prec, int Size>
    [[nodiscard]] __device__ auto reciprocal(device_obj<SIMD<Real<Prec>, Size>> x) noexcept {
        using T = Real<Prec>;
        using ResultType = device_obj<SIMD<T, Size>>;
        return ResultType(T(1)) / x;
    }
}
