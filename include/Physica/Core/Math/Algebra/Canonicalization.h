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
#include "Physica/Core/Scalar/ExprID.h"
#include "Physica/Core/Scalar/Scalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/Tensor.h"

namespace Physica {
    template<ExprID, class LHS, class RHS = LHS> class VectorExpr;
    template<ExprID, class LHS, class RHS = LHS> class MatrixExpr;
    template<ExprID, class LHS, class RHS = LHS> class TensorExpr;

    namespace Internal {
        template<class T>
        consteval int concept_score() noexcept {
            // Increasing order; detailed number is not important
            if (Scalar<T>)
                return 1;
            if (Packet<T>)
                return 2;
            if (Vector<T>)
                return 3;
            if (Matrix<T>)
                return 4;
            if (Tensor<T>)
                return 5;
            return 0;
        }

        template<class T>
        consteval int64_t canonical_score() noexcept {
            int64_t result = 0;
            if constexpr (!Scalar<T> && !Packet<T>)
                result -= int64_t(T::isSparse());
            result += int64_t(T::isComplex());
            result += int64_t(T::isDiffable()) * 2;
            result += int64_t(instanceof_xt<VectorExpr, T>) * 10;
            result += int64_t(instanceof_xt<MatrixExpr, T>) * 10;
            result += int64_t(instanceof_xt<TensorExpr, T>) * 10;
            return result;
        }

        template<class T>
        consteval bool canonicalizable() noexcept {
            return concept_score<T>() > 0;
        }

        consteval bool canonicalizable(const auto& expr) noexcept {
            return canonicalizable<decltype(expr)>();
        }
    }

    template<class T1, class T2>
    consteval bool canonicalized() {
        using namespace Internal;
        using U1 = std::remove_cvref_t<T1>;
        using U2 = std::remove_cvref_t<T2>;
        if (std::same_as<U1, U2>)
            return true;
        if (concept_score<T1>() == concept_score<T2>())
            return canonical_score<U1>() >= canonical_score<std::remove_cvref_t<U2>>(); // TODO: Maybe give a warning if scores equal
        return concept_score<T1>() > concept_score<T2>();
    }

    consteval bool canonicalized(const auto& expr1, const auto& expr2) noexcept {
        using namespace Internal;
        static_assert(canonicalizable(expr1) && canonicalizable(expr2));
        return canonicalized<decltype(expr1), decltype(expr2)>();
    }
}
