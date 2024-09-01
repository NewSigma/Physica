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
    class VectorExpr<ExpressionType::Exp, VectorType> : public UnitaryVectorExpr<ExpressionType::Exp, VectorType> {
        using This = VectorExpr<ExpressionType::Exp, VectorType>;
        using Base = UnitaryVectorExpr<ExpressionType::Exp, VectorType>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const { return exp(Base::getExpr().calc(index)); }
    };

    template<class VectorType>
    [[nodiscard]] inline auto exp(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Exp, VectorType>(v);
    }
}
