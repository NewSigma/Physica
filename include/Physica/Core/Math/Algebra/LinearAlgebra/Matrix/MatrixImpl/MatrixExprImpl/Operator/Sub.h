/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
