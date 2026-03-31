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
    class VectorExpr<ExprID::Ln, V> : public UnitaryVectorExpr<ExprID::Ln, V> {
        using This = VectorExpr<ExprID::Ln, V>;
        using Base = UnitaryVectorExpr<ExprID::Ln, V>;
    public:
        using Base::isReverseDiff;
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
        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        void reverse(const auto& grad) const noexcept;
        [[nodiscard]] auto values(this auto&&) noexcept;
};

    template<Vector V>
    auto VectorExpr<ExprID::Ln, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return ln(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Ln, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        auto x = input.template load<Size>();
        assert(x.isPositive().horizontal_and());
        return ln(x);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Ln, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return ln(input.template load<Size>(count)).cutoff(count);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Ln, V>::calc(size_t index) const -> CoDiff<T> {
        return ln(Base::getExpr().calc(index));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Ln, V>::calc_value(size_t index) const -> Tv {
        return ln(Base::getExpr().calc_value(index));
    }

    template<Vector V>
    void VectorExpr<ExprID::Ln, V>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& expr = Base::getExpr();
        if constexpr (Scalar<decltype(grad)>)
            expr.reverse(grad.value() / expr.values());
        else {
            static_assert(Vector<decltype(grad)>, "[Error]: Unexpected type");
            expr.reverse(divide(grad.values(), expr.values()));
        }
    }

    template<Vector V>
    auto VectorExpr<ExprID::Ln, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return ln(std::forward<Self>(self).getExpr().values());
    }

template<Vector V>
    [[nodiscard, gnu::always_inline]] auto ln(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Ln, V&&>(std::forward<V>(v));
    }
}
