/*
 * Copyright 2025 Weibo He.
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

#include "../MatrixExpr.cuh"

namespace Physica {
    template<Matrix T>
    class device_obj<MatrixExpr<ExprType::Sqrt, T>>
            : public device_obj<UnitaryMatrixExpr<ExprType::Sqrt, T>> {
        using Base = device_obj<UnitaryMatrixExpr<ExprType::Sqrt, T>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return sqrt(Base::getExpr().calc(row, col)); }

        [[nodiscard]] __device__ ValueType calc_value(size_t row, size_t col) const {
            return sqrt(Base::getExpr().calc_value(row, col));
        }
    };

    template<Matrix T>
    [[nodiscard]] __host__ __device__ inline auto sqrt_elem(T&& m) noexcept requires(CUDA<T>) {
        return device_obj<MatrixExpr<ExprType::Sqrt, T&&>>(std::forward<T>(m));
    }
}
