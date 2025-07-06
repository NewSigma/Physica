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
    class VectorExpr<ExprType::Tanh, V> : public UnitaryVectorExpr<ExprType::Tanh, V> {
        using This = VectorExpr<ExprType::Tanh, V>;
        using Base = UnitaryVectorExpr<ExprType::Tanh, V>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return tanh(Base::getExpr().calc(index)); }

        [[nodiscard]] Tv calc_value(size_t index) const { return tanh(Base::getExpr().calc_value(index)); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return tanh(Base::getExpr().template packet<Pack>(index));
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return tanh(Base::getExpr().template packetPartial<Pack>(index, count));
        }

        void reverse(const Vector auto& y, const Vector auto& grad_) const noexcept requires(isReverseDiff);
    };

    template<Vector V>
    void VectorExpr<ExprType::Tanh, V>::reverse(const Vector auto& y, const Vector auto& grad_) const noexcept requires(isReverseDiff) {
        const auto& expr = Base::getExpr();
        expr.reverse(hadamard((Tv(1) - square(y.values())), grad_));
    }

    template<Vector V>
    [[nodiscard]] inline auto tanh(V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Tanh, V&&>(std::forward<V>(v));
    }
}
