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
    class VectorExpr<ExprID::Tanh, V> : public UnitaryVectorExpr<ExprID::Tanh, V> {
        using This = VectorExpr<ExprID::Tanh, V>;
        using Base = UnitaryVectorExpr<ExprID::Tanh, V>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        template<Packet Pack>
        [[nodiscard]] static Pack operator()(std::random_access_iterator auto input) noexcept;
        template<Packet Pack>
        [[nodiscard]] static Pack operator()(std::random_access_iterator auto input, size_t count) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return tanh(Base::getExpr().calc(index)); }
        [[nodiscard]] Tv calc_value(size_t index) const { return tanh(Base::getExpr().calc_value(index)); }

        void reverse(const Vector auto& y, const Vector auto& grad) const noexcept;
    };

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Tanh, V>::operator()(std::random_access_iterator auto input) noexcept {
        return tanh(input.template load<Pack>());
    }

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Tanh, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept {
        return tanh(input.template load<Pack>(count));
    }

    template<Vector V>
    void VectorExpr<ExprID::Tanh, V>::reverse(const Vector auto& y, const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& expr = Base::getExpr();
        expr.reverse(hadamard((Tv(1) - square(y.values())), grad));
    }

    template<Vector V>
    [[nodiscard]] auto tanh(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Tanh, V&&>(std::forward<V>(v));
    }
}
