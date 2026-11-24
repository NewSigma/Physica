/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
