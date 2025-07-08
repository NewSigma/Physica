/*
 * Copyright 2024-2025 Weibo He.
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
    class VectorExpr<ExprType::Relu, V> : public UnitaryVectorExpr<ExprType::Relu, V> {
        using This = VectorExpr<ExprType::Relu, V>;
        using Base = UnitaryVectorExpr<ExprType::Relu, V>;
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
        void reverse(const Vector auto& grad_) const noexcept requires(isReverseDiff);
    };

    template<Vector V>
    void VectorExpr<ExprType::Relu, V>::reverse(const Vector auto& grad_) const noexcept requires(isReverseDiff) {
        const auto& grad = grad_.values();
        for (size_t i = 0; i < Base::getLength(); ++i) {
            auto v = Base::getExpr().calc(i);
            v.reverse(v.isPositive() ? grad.calc(i) : Tv(0));
        }
    }

    template<Vector V>
    [[nodiscard]] auto relu(V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Relu, V&&>(std::forward<V>(v));
    }
}
