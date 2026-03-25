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
    template<Vector V, Scalar U>
    class VectorExpr<ExprID::Pow, V, U>
            : public BinaryVectorExpr<ExprID::Pow, V, U> {
        using Base = BinaryVectorExpr<ExprID::Pow, V, U>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs, size_t count) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t i) const;
        [[nodiscard]] Tv calc_value(size_t i) const;

        void reverse(const Vector auto& grad) const noexcept;
    };

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Pow, V, U>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept -> CoDiff<T> {
        return pow(*lhs, rhs);
    }

    template<Vector V, Scalar U>
    template<int Size>
    auto VectorExpr<ExprID::Pow, V, U>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept -> SIMD<T, Size> {
        return pow(lhs.template load<Size>(), SIMD<T, Size>(rhs));
    }

    template<Vector V, Scalar U>
    template<int Size>
    auto VectorExpr<ExprID::Pow, V, U>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs, size_t count) noexcept -> SIMD<T, Size> {
        return pow(lhs.template load<Size>(count), SIMD<T, Size>(rhs));
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Pow, V, U>::calc(size_t i) const -> CoDiff<T> {
        return pow(Base::getLHS().calc(i), Base::getRHS());
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Pow, V, U>::calc_value(size_t i) const -> Tv {
        return pow(Base::getLHS().calc_value(i), Base::getRHS().value());
    }

    template<Vector V, Scalar U>
    [[nodiscard, gnu::always_inline]] auto pow(V&& v, U&& x) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Pow, V&&, U&&>(std::forward<V>(v), std::forward<U>(x));
    }
}
