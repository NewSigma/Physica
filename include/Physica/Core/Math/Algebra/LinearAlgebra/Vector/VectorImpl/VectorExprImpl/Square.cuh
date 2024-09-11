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
    class device_obj<VectorExpr<ExprType::Square, VectorType>>
            : public device_obj<UnitaryVectorExpr<ExprType::Square, VectorType>> {
        using Base = device_obj<UnitaryVectorExpr<ExprType::Square, VectorType>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t s) const {
            return square(Base::template getExpr<Owner>().template calc<Owner>(s));
        }
    };

    template<class VectorType>
    [[nodiscard]] __host__ __device__ inline auto square(const device_obj<RValueVector<VectorType>>& v) noexcept {
        return device_obj<VectorExpr<ExprType::Square, VectorType>>(v);
    }
}
