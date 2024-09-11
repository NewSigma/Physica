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
    class device_obj<VectorExpr<ExprType::Div, VectorType, ScalarBase<AnyScalar>>>
            : public device_obj<BinaryVectorExpr<ExprType::Div, VectorType, ScalarBase<AnyScalar>>> {
        using Base = device_obj<BinaryVectorExpr<ExprType::Div, VectorType, ScalarBase<AnyScalar>>>;
        static_assert(is_scalar<AnyScalar>::value, "[Error]: This is not a scalar type");
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return ScalarType(Base::template getLHS<Owner>().template calc<Owner>(index)) / ScalarType(Base::template getRHS<Owner>());
        }
    };

    template<class VectorType, class ScalarType>
    [[nodiscard]] __host__ __device__ inline device_obj<VectorExpr<ExprType::Div, VectorType, ScalarBase<ScalarType>>>
    operator/(const device_obj<RValueVector<VectorType>>& v, const ScalarBase<ScalarType>& s) noexcept {
        return {v.getDerived(), s.getDerived()};
    }
}
