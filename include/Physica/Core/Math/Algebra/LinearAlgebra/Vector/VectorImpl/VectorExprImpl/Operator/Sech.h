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
    class VectorExpr<ExprID::Sech, V> : public UnitaryVectorExpr<ExprID::Sech, V> {
        using This = VectorExpr<ExprID::Sech, V>;
        using Base = UnitaryVectorExpr<ExprID::Sech, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto input) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return sech(Base::getExpr().calc(index)); }

        [[nodiscard]] Tv calc_value(size_t index) const { return sech(Base::getExpr().calc_value(index)); }
        [[nodiscard]] auto values(this auto&&) noexcept;
};

    template<Vector V>
    auto VectorExpr<ExprID::Sech, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return sech(*input);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Sech, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return sech(std::forward<Self>(self).getExpr().values());
    }

template<Vector V>
    [[nodiscard, gnu::always_inline]] auto sech(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Sech, V&&>(std::forward<V>(v));
    }
}
