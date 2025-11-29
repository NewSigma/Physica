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
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;

        void reverse(const auto& grad) const noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Ln, V>::calc(size_t index) const -> CoDiff<T> {
        return ln(Base::getExpr().calc(index));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Ln, V>::calc_value(size_t index) const -> Tv {
        return ln(Base::getExpr().calc_value(index));
    }

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Ln, V>::packet(size_t index) const {
        auto x = Base::getExpr().template packet<Pack>(index);
        assert(x.isPositive().horizontal_and());
        return ln(x);
    }

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Ln, V>::packetPartial(size_t index, size_t count) const {
        return ln(Base::getExpr().template packetPartial<Pack>(index, count)).cutoff(count);
    }

    template<Vector V>
    void VectorExpr<ExprID::Ln, V>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& expr = Base::getExpr();
        if constexpr (Scalar<decltype(grad)>)
            expr.reverse(grad.value() / expr.values());
        else {
            static_assert(Vector<decltype(grad)>, "[Error]: Unexpected type");
            expr.reverse(divide(grad.values(), expr.values()));
        }
    }

    template<Vector V>
    [[nodiscard]] auto ln(V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprID::Ln, V&&>(std::forward<V>(v));
    }
}
