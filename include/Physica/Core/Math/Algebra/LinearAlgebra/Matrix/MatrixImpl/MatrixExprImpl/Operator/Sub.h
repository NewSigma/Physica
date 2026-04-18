/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.h"

namespace Physica {
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator-(Matrix auto&& m, Scalar auto&& x) noexcept {
        return std::forward<decltype(m)>(m) + (-std::forward<decltype(x)>(x));
    }

    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator-(Scalar auto&& x, Matrix auto&& m) noexcept {
        return std::forward<decltype(x)>(x) + (-std::forward<decltype(m)>(m));
    }

    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator-(Matrix auto&& m1, Matrix auto&& m2) noexcept {
        return std::forward<decltype(m1)>(m1) + (-std::forward<decltype(m2)>(m2));
    }
}
