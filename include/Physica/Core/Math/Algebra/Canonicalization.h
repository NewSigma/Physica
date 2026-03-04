/*
 * Copyright 2026 Weibo He.
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

#include <type_traits>
#include "Physica/Core/Scalar/Scalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/Tensor.h"

namespace Physica {
    consteval bool canonicalizable(const auto& expr) noexcept {
        using E = decltype(expr);
        return Scalar<E> || Vector<E> || Matrix<E> || Tensor<E>;
    }

    consteval bool canonicalized(const auto& expr1, const auto& expr2) noexcept {
        static_assert(canonicalizable(expr1) && canonicalizable(expr2));

        using E1 = std::remove_cvref_t<decltype(expr1)>;
        using E2 = std::remove_cvref_t<decltype(expr2)>;
        if constexpr (!E1::isComplex && E2::isComplex)
            return false;
        if constexpr (!E1::isDiffable && E2::isDiffable)
            return false;
        return true;
    }
}
