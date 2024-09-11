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

namespace Physica::Core {
    template<class VectorType>
    class VectorExpr<ExprType::Cos, VectorType> : public UnitaryVectorExpr<ExprType::Cos, VectorType> {
        using This = VectorExpr<ExprType::Cos, VectorType>;
        using Base = UnitaryVectorExpr<ExprType::Cos, VectorType>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const { return cos(Base::getExpr().calc(index)); }
    };

    template<class VectorType>
    [[nodiscard]] inline auto cos(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExprType::Cos, VectorType>(v);
    }
}
