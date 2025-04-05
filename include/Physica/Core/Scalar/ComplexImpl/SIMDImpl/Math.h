/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Scalar/ComplexImpl/SIMD.h"

namespace Physica {
    template<Scalar T, size_t Size>
    SIMD<T, Size * 2> abs(const SIMD<Complex<T>, Size>& x);

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size> sqrt(const SIMD<Complex<T>, Size>& x);

    template<Scalar T, size_t Size>
    SIMD<Complex<T>, Size> exp(const SIMD<Complex<T>, Size>& x);

    template<Scalar T, size_t Size>
    SIMD<Complex<T>, Size> ln(const SIMD<Complex<T>, Size>& x);

    template<Scalar T, size_t Size>
    SIMD<Complex<T>, Size> ln1p(const SIMD<Complex<T>, Size>& x);

    template<Scalar T, size_t Size>
    SIMD<Complex<T>, Size> tanh(const SIMD<Complex<T>, Size>& x);

    template<Scalar T, size_t Size>
    SIMD<Complex<T>, Size> lncosh(const SIMD<Complex<T>, Size>& x);
}

#include "MathImpl.h"
