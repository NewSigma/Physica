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
    class VectorExpr<ExprID::Ln1p, V> : public UnitaryVectorExpr<ExprID::Ln1p, V> {
        using This = VectorExpr<ExprID::Ln1p, V>;
        using Base = UnitaryVectorExpr<ExprID::Ln1p, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Trv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto input) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input, size_t count) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return ln1p(Base::getExpr().calc(index)); }
        [[nodiscard]] Tv calc_value(size_t index) const { return ln1p(Base::getExpr().calc_value(index)); }

        void reverse(const Scalar auto& grad) const noexcept;
        [[nodiscard]] auto values(this auto&&) noexcept;
};

    template<Vector V>
    auto VectorExpr<ExprID::Ln1p, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return ln1p(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Ln1p, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return ln1p(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Ln1p, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return ln1p(input.template load<Size>(count));
    }

    template<Vector V>
    void VectorExpr<ExprID::Ln1p, V>::reverse(const Scalar auto& grad) const noexcept {
        static_assert(Base::isReverseDiff);
        const auto& expr = Base::getExpr();
        expr.reverse(divide(grad.value(), (Trv(1) + expr.values())));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Ln1p, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return ln1p(std::forward<Self>(self).getExpr().values());
    }

template<Vector V>
    [[nodiscard, gnu::always_inline]] auto ln1p(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Ln1p, V&&>(std::forward<V>(v));
    }
}
