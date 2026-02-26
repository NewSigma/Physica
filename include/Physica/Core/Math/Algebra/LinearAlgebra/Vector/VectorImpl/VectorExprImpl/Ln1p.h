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
    class VectorExpr<ExprID::Ln1p, V> : public UnitaryVectorExpr<ExprID::Ln1p, V> {
        using This = VectorExpr<ExprID::Ln1p, V>;
        using Base = UnitaryVectorExpr<ExprID::Ln1p, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Trv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return ln1p(Base::getExpr().calc(index)); }
        [[nodiscard]] Tv calc_value(size_t index) const { return ln1p(Base::getExpr().calc_value(index)); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const noexcept;
        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index, size_t count) const noexcept;

        void reverse(const Scalar auto& grad) const noexcept;
    };

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Ln1p, V>::packet(size_t index) const noexcept {
        return ln1p(Base::getExpr().template packet<Pack>(index));
    }

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Ln1p, V>::packet(size_t index, size_t count) const noexcept {
        return ln1p(Base::getExpr().template packet<Pack>(index, count));
    }

    template<Vector V>
    void VectorExpr<ExprID::Ln1p, V>::reverse(const Scalar auto& grad) const noexcept {
        static_assert(Base::isReverseDiff);
        const auto& expr = Base::getExpr();
        expr.reverse(divide(grad.value(), (Trv(1) + expr.values())));
    }

    template<Vector V>
    [[nodiscard]] auto ln1p(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Ln1p, V&&>(std::forward<V>(v));
    }
}
