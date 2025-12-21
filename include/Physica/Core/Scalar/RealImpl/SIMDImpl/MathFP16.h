/*
 * Copyright 2025 Weibo He.
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
