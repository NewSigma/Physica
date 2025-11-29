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
    template<Vector T>
    class device_obj<VectorExpr<ExprID::Exp, T>>
            : public device_obj<UnitaryVectorExpr<ExprID::Exp, T>> {
        using Base = device_obj<UnitaryVectorExpr<ExprID::Exp, T>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return exp(Base::getExpr().calc(index));
        }
    };

    template<Vector T>
    [[nodiscard]] __host__ __device__ auto exp(T&& v) noexcept requires(CUDA<T>) {
        return device_obj<VectorExpr<ExprID::Exp, T&&>>(std::forward<T>(v));
    }
}
