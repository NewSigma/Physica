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
    template<Matrix, Vector> class GEMV;

    template<Vector V>
    class VectorExpr<ExprID::Minus, V>
            : public UnitaryVectorExpr<ExprID::Minus, V> {
        using Base = UnitaryVectorExpr<ExprID::Minus, V>;
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
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const;

        [[nodiscard]] CoDiff<T> calc(size_t s) const;

        void reverse(const auto& grad) const noexcept;
        [[nodiscard]] auto values(this auto&&) noexcept;
};

    template<Vector V>
    auto VectorExpr<ExprID::Minus, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return -(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Minus, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return -input.template load<Size>();
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Minus, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return -input.template load<Size>(count);
    }

    template<Vector V>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Minus, V>::assign(Vector auto&& v) const {
        const auto& expr = Base::getExpr();
        constexpr bool IsBinaryExpr = requires { expr.getLHS() * expr.getRHS(); };

        if constexpr (IsBinaryExpr) {
            const auto& lhs = expr.getLHS();
            const auto& rhs = expr.getRHS();
            constexpr bool IsMatVecProd = Matrix<decltype(lhs)> && Vector<decltype(rhs)>;
            if constexpr (IsMatVecProd) {
                if constexpr (Traits<std::remove_cvref_t<decltype(lhs)>>::FastAssign)
                    ((-lhs) * rhs).template assign<P>(v);
                else
                    (lhs * (-rhs)).template assign<P>(v);
            }
            else
                Base::template assign<P>(v);
        }
        else
            Base::template assign<P>(v);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Minus, V>::calc(size_t s) const -> CoDiff<T> {
        return -Base::getExpr().calc(s);
    }

    template<Vector V>
    void VectorExpr<ExprID::Minus, V>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff());
        Base::getExpr().reverse(-grad);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Minus, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return -std::forward<Self>(self).getExpr().values();
    }

template<Vector V>
    [[nodiscard, gnu::always_inline]] auto operator-(V&& v) noexcept requires(!DeviceObj<V>) {
        if constexpr (instanceof<GEMV, V>)
            return v.getLHS() * (-v.getRHS());
        else
            return VectorExpr<ExprID::Minus, V&&>(std::forward<V>(v));
    }
}
