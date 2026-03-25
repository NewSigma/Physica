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
    class VectorExpr<ExprID::Sin, V> : public UnitaryVectorExpr<ExprID::Sin, V> {
        using This = VectorExpr<ExprID::Sin, V>;
        using Base = UnitaryVectorExpr<ExprID::Sin, V>;
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
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return sin(Base::getExpr().calc(index)); }
        [[nodiscard]] Tv calc_value(size_t index) const { return sin(Base::getExpr().calc_value(index)); }
    };

    template<Vector V>
    auto VectorExpr<ExprID::Sin, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return sin(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Sin, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return sin(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Sin, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return sin(input.template load<Size>(count));
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto sin(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Sin, V&&>(std::forward<V>(v));
    }
}
