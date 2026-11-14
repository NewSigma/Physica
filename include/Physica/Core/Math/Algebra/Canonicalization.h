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

#include <compare>
#include <string_view>
#include <type_traits>
#include "Physica/Core/Math/Algebra/ExprID.h"
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
            using U = std::remove_cvref_t<remove_device_obj_t<std::remove_cvref_t<T>>>;
            // Increasing order; detailed number is not important
            if (Scalar<U>)
                return 1;
            if (Packet<U>)
                return 2;
            if (Vector<U>)
                return 3;
            if (Matrix<U>)
                return 4;
            if (Tensor<U>)
                return 5;
            return 0;
        }

        template<class T>
        consteval int64_t canonical_score() noexcept {
            using U = std::remove_cvref_t<remove_device_obj_t<std::remove_cvref_t<T>>>;
            int64_t result = 0;
            if constexpr (!Scalar<U> && !Packet<U>)
                result -= int64_t(U::isSparse());
            result += int64_t(U::isComplex());
            result += int64_t(U::isDiffable()) * 2;
            result += int64_t(instanceof_x<U, VectorExpr>) * 10;
            result += int64_t(instanceof_x<U, MatrixExpr>) * 10;
            result += int64_t(instanceof_x<U, TensorExpr>) * 10;
            return result;
        }

        template<class T>
        consteval bool canonicalizable() noexcept {
            return concept_score<T>() > 0;
        }

        consteval bool canonicalizable(const auto& expr) noexcept {
            return canonicalizable<decltype(expr)>();
        }

        struct CanonicalKey {
            int conceptScore = 0;
            int64_t score = 0;

            constexpr auto operator<=>(const CanonicalKey& other) const noexcept = default;
        };

        template<class T>
        consteval auto canonical_key() noexcept {
            using U = std::remove_cvref_t<T>;
            return CanonicalKey{concept_score<U>(), canonical_score<U>()};
        }

        template<class T>
        consteval std::string_view pretty_typename() noexcept {
            constexpr std::string_view marker = "T = ";
            const std::string_view pretty = __PRETTY_FUNCTION__;
            const auto pos = pretty.find(marker);
            bool isUnrecognized = pos == std::string_view::npos;
            if (isUnrecognized)
                return pretty; // Still unique per T

            std::string_view result = pretty.substr(pos + marker.size());
            if (!result.empty() && result.back() == ']')
                result.remove_suffix(1);
            return result;
        }

        template<class T, class U>
        consteval std::strong_ordering pretty_type_order() noexcept {
            const int cmp = pretty_typename<T>().compare(pretty_typename<U>());
            if (cmp < 0)
                return std::strong_ordering::less;
            if (cmp > 0)
                return std::strong_ordering::greater;
            return std::strong_ordering::equal;
        }

        template<class T, class U>
        consteval std::strong_ordering type_order() noexcept {
        #if defined(__cpp_lib_type_order) // CXX26
            return std::type_order_v<T, U>;
        #elifndef __has_builtin
            return pretty_type_order<T, U>();
        #elif __has_builtin(__builtin_type_order)
            return __builtin_type_order(T, U);
        #else
            return pretty_type_order<T, U>(); // TODO: Drop this branch clang>=24 or gcc>=16 is required
        #endif
        }
    }

    template<class T1, class T2>
    consteval bool canonicalized() {
        using namespace Internal;
        using U1 = std::remove_cvref_t<T1>;
        using U2 = std::remove_cvref_t<T2>;
        if (std::same_as<U1, U2>)
            return true;
        const auto key1 = canonical_key<U1>();
        const auto key2 = canonical_key<U2>();
        if (key1 != key2)
            return key1 > key2;
        return type_order<U1, U2>() > 0;
    }

    consteval bool canonicalized(const auto& expr1, const auto& expr2) noexcept {
        using namespace Internal;
        static_assert(canonicalizable(expr1) && canonicalizable(expr2));
        return canonicalized<decltype(expr1), decltype(expr2)>();
    }
}
