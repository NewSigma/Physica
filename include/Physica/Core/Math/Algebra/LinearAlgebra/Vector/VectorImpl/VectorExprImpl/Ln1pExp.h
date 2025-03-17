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
    class VectorExpr<ExprType::Ln1pExp, V> : public UnitaryVectorExpr<ExprType::Ln1pExp, V> {
        using This = VectorExpr<ExprType::Ln1pExp, V>;
        using Base = UnitaryVectorExpr<ExprType::Ln1pExp, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t index) const { return ln1pexp(Base::getExpr().calc(index)); }

        [[nodiscard]] Tv calc_value(size_t index) const { return ln1pexp(Base::getExpr().calc_value(index)); }
    };

    template<Vector V>
    [[nodiscard]] inline auto ln1pexp(V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Ln1pExp, V&&>(std::forward<V>(v));
    }
}
