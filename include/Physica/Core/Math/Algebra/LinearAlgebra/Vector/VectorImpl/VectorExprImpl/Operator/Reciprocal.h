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
    class VectorExpr<ExprID::Reciprocal, V> : public UnitaryVectorExpr<ExprID::Reciprocal, V> {
        using This = VectorExpr<ExprID::Reciprocal, V>;
        using Base = UnitaryVectorExpr<ExprID::Reciprocal, V>;
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
        [[nodiscard]] Pack packet(size_t index) const noexcept;
        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index, size_t count) const noexcept;

        void reverse(const Vector auto& grad) const noexcept;
        void reverse(const Vector auto& y, const Vector auto& grad) const noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Reciprocal, V>::calc(size_t index) const -> CoDiff<T> {
        return reciprocal(Base::getExpr().calc(index));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Reciprocal, V>::calc_value(size_t index) const -> Tv { return reciprocal(Base::getExpr().calc_value(index)); }

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Reciprocal, V>::packet(size_t index) const noexcept {
        auto x = Base::getExpr().template packet<Pack>(index);
        assert(!x.isZero().horizontal_or() && "[Error]: Divide by zero");
        return reciprocal(x);
    }

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Reciprocal, V>::packet(size_t index, size_t count) const noexcept {
        auto x = reciprocal(Base::getExpr().template packet<Pack>(index, count)).cutoff(count);
        assert(x.isFinite().horizontal_and() && "[Error]: Divide by zero");
        return x;
    }

    template<Vector V>
    void VectorExpr<ExprID::Reciprocal, V>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        reverse(Base::getExpr().values(), grad);
    }

    template<Vector V>
    void VectorExpr<ExprID::Reciprocal, V>::reverse(const Vector auto& y, const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        Base::getExpr().reverse(hadamard(-square(y.values()), grad));
    }

    template<Vector V>
    [[nodiscard]] auto reciprocal(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Reciprocal, V&&>(std::forward<V>(v));
    }
}
