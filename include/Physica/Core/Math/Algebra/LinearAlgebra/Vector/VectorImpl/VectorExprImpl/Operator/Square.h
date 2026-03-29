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
    class VectorExpr<ExprID::Square, V> : public UnitaryVectorExpr<ExprID::Square, V> {
        using This = VectorExpr<ExprID::Square, V>;
        using Base = UnitaryVectorExpr<ExprID::Square, V>;
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
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const;
        void assign_mkl(Vector auto& v) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        void reverse(const auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        using Base::getExpr;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Square, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return square(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Square, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return square(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Square, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return square(input.template load<Size>(count));
    }

    template<Vector V>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Square, V>::assign(Vector auto&& v) const {
        Base::template assign_base<P>(v);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Square, V>::calc(size_t index) const -> CoDiff<T> {
        return square(getExpr().calc(index));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Square, V>::calc_value(size_t index) const -> Tv {
        return square(getExpr().calc_value(index));
    }

    template<Vector V>
    void VectorExpr<ExprID::Square, V>::reverse(const auto& grad) const noexcept {
        static_assert(Base::isReverseDiff);
        using G = decltype(grad);
        const auto& expr = getExpr();
        if constexpr (Scalar<G>)
            expr.reverse(expr.values() * (Trv(2) * grad));
        else {
            static_assert(Vector<G>, "[Error]: Unexpected type");
            expr.reverse(hadamard(expr.values(), grad) * Trv(2));
        }
    }

    template<Vector V>
    auto VectorExpr<ExprID::Square, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return square(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto square(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Square, V&&>(std::forward<V>(v));
    }
}

#ifdef PHYSICA_MKL
    #include "MKL/Square.h"
#endif
