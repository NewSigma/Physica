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

#include "../VectorExpr.h"

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
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const {
            return relu(Base::getExpr().calc(index));
        }

        [[nodiscard]] Tv calc_value(size_t index) const {
            return relu(Base::getExpr().calc_value(index));
        }

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;
    };

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
