/*
 * Copyright 2024 Weibo He.
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

namespace Physica {
    template<Vector T>
    class VectorExpr<ExprType::Sech, T> : public UnitaryVectorExpr<ExprType::Sech, T> {
        using This = VectorExpr<ExprType::Sech, T>;
        using Base = UnitaryVectorExpr<ExprType::Sech, T>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t index) const { return sech(Base::getExpr().calc(index)); }

        [[nodiscard]] ValueType calc_value(size_t index) const { return sech(Base::getExpr().calc_value(index)); }
    };

    template<Vector T>
    [[nodiscard]] inline auto sech(const T& v) noexcept {
        return VectorExpr<ExprType::Sech, T>(v);
    }
}
