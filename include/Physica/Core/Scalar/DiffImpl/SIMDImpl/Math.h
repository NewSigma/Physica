/*
 * Copyright 2024 Weibo He.
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

#include "Physica/Core/Scalar/Scalar.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto abs(const SIMD<Diff<T, Mode, Order>, Size>& x);

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto square(const SIMD<Diff<T, Mode, Order>, Size>& x);

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto reciprocal(const SIMD<Diff<T, Mode, Order>, Size>& x);

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto ln(const SIMD<Diff<T, Mode, Order>, Size>& x);

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto ln1p(const SIMD<Diff<T, Mode, Order>, Size>& x);

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto exp(const SIMD<Diff<T, Mode, Order>, Size>& x);

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto sincos(const SIMD<Diff<T, Mode, Order>, Size>& x, SIMD<Diff<T, Mode, Order>, Size>& s, SIMD<Diff<T, Mode, Order>, Size>& c);

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto tanh(const SIMD<Diff<T, Mode, Order>, Size>& x);

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto lncosh(const SIMD<Diff<T, Mode, Order>, Size>& x);
}

#include "MathImpl.h"
