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
    class VectorExpr<ExprID::Relu, V> : public UnitaryVectorExpr<ExprID::Relu, V> {
        using This = VectorExpr<ExprID::Relu, V>;
        using Base = UnitaryVectorExpr<ExprID::Relu, V>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input, size_t count) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;
    };

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Relu, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return relu(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Relu, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return relu(input.template load<Size>(count));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Relu, V>::calc(size_t index) const -> CoDiff<T> {
        return relu(Base::getExpr().calc(index));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Relu, V>::calc_value(size_t index) const -> Tv {
        return relu(Base::getExpr().calc_value(index));
    }

    template<Vector V>
    void VectorExpr<ExprID::Relu, V>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& g = grad.values();
        for (size_t i = 0; i < Base::getLength(); ++i) {
            auto v = Base::getExpr().calc(i);
            v.reverse(v.isPositive() ? g.calc(i) : Tv(0));
        }
    }

    template<Vector V>
    [[nodiscard]] auto relu(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Relu, V&&>(std::forward<V>(v));
    }
}
