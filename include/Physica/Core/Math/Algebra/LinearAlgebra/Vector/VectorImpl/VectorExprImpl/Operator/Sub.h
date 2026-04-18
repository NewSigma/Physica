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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.h"

namespace Physica {
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator-(Vector auto&& v, Scalar auto&& x) noexcept {
        return std::forward<decltype(v)>(v) + (-std::forward<decltype(x)>(x));
    }

    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator-(Scalar auto&& x, Vector auto&& v) noexcept {
        return std::forward<decltype(x)>(x) + (-std::forward<decltype(v)>(v));
    }

    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator-(Vector auto&& v1, Vector auto&& v2) noexcept {
        return std::forward<decltype(v1)>(v1) + (-std::forward<decltype(v2)>(v2));
    }
}
