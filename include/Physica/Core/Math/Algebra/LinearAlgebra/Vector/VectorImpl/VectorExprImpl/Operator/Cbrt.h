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
    template<Vector V>
    class VectorExpr<ExprID::Cbrt, V> : public UnitaryVectorExpr<ExprID::Cbrt, V> {
        using This = VectorExpr<ExprID::Cbrt, V>;
        using Base = UnitaryVectorExpr<ExprID::Cbrt, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto input) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input, size_t count) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t s) const { return cbrt(Base::getExpr().calc(s)); }
        [[nodiscard]] Tv calc_value(size_t index) const { return cbrt(Base::getExpr().calc_value(index)); }
    };

    template<Vector V>
    auto VectorExpr<ExprID::Cbrt, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return cbrt(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Cbrt, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        SIMD<T, Size> result = input.template load<Size>();
        for (int i = 0; i < Size; ++i)
            result.insert(i, cbrt(result[i]));
        return result;
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Cbrt, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        SIMD<T, Size> result = input.template packet<Size>(count);
        for (size_t i = 0; i < count; ++i)
            result.insert(i, cbrt(result[i]));
        return result;
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto cbrt(V&& v) noexcept {
        return VectorExpr<ExprID::Cbrt, V&&>(std::forward<V>(v));
    }
}
