/*
 * Copyright 2025-2026 Weibo He.
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
    class VectorExpr<ExprID::ArcTanh, V> : public UnitaryVectorExpr<ExprID::ArcTanh, V> {
        using This = VectorExpr<ExprID::ArcTanh, V>;
        using Base = UnitaryVectorExpr<ExprID::ArcTanh, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto input) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return arctanh(Base::getExpr().calc(index)); }

        [[nodiscard]] auto values(this auto&&) noexcept;
};

    template<Vector V>
    auto VectorExpr<ExprID::ArcTanh, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return arctanh(*input);
    }

    template<Vector V>
    auto VectorExpr<ExprID::ArcTanh, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return arctanh(std::forward<Self>(self).getExpr().values());
    }

template<Vector V>
    [[nodiscard, gnu::always_inline]] auto arctanh(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::ArcTanh, V&&>(std::forward<V>(v));
    }
}
