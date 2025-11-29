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

#include "../VectorExpr.cuh"

namespace Physica {
    template<Vector V>
    class device_obj<VectorExpr<ExprID::Minus, V>>
            : public device_obj<UnitaryVectorExpr<ExprID::Minus, V>> {
        using Base = device_obj<UnitaryVectorExpr<ExprID::Minus, V>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t index) const;
        [[nodiscard]] __device__ Tv calc_value(size_t index) const;
    };

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprID::Minus, V>>::calc(size_t index) const -> T {
        if constexpr (isReverseDiff)
            return calc_value(index);
        else
            return -Base::getExpr().calc(index);
    }

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprID::Minus, V>>::calc_value(size_t index) const -> Tv {
        return -Base::getExpr().calc_value(index);
    }

    template<Vector V>
    [[nodiscard]] __host__ __device__ auto operator-(V&& v) noexcept requires(CUDA<V>) {
        return device_obj<VectorExpr<ExprID::Minus, V&&>>(std::forward<V>(v));
    }
}
