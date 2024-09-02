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
    template<class VectorType, class AnyScalar>
    class VectorExpr<ExpressionType::Pow, VectorType, ScalarBase<AnyScalar>>
            : public BinaryVectorExpr<ExpressionType::Pow, VectorType, ScalarBase<AnyScalar>> {
        using Base = BinaryVectorExpr<ExpressionType::Pow, VectorType, ScalarBase<AnyScalar>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return pow(Base::getLHS().calc(i), Base::getRHS()); }
    };

    template<class VectorType, class AnyScalar>
    [[nodiscard]] inline auto pow(const RValueVector<VectorType>& v, const ScalarBase<AnyScalar>& s) noexcept {
        return VectorExpr<ExpressionType::Pow, VectorType, ScalarBase<AnyScalar>>(v.getDerived(), s.getDerived());
    }
}
