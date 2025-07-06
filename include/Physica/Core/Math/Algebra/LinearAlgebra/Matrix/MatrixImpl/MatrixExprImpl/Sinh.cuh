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
    template<Matrix M>
    class device_obj<MatrixExpr<ExprType::Sinh, M>>
            : public device_obj<UnitaryMatrixExpr<ExprType::Sinh, M>> {
        using Base = device_obj<UnitaryMatrixExpr<ExprType::Sinh, M>>;
    public:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const {
            if constexpr (isReverseDiff)
                return calc_value(row, col);
            else
                return sinh(Base::getExpr().calc(row, col));
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            return sinh(Base::getExpr().calc_value(row, col));
        }

        void reverse(const Matrix auto& grad) const noexcept requires(isReverseDiff);
        using Base::reverse;
    };

    template<Matrix M>
    void device_obj<MatrixExpr<ExprType::Sinh, M>>::reverse(const Matrix auto& grad) const noexcept requires(isReverseDiff) {
        const auto& expr = Base::getExpr();
        expr.reverse(hadamard(cosh_elem(expr.values()), grad));
    }

    template<Matrix M>
    [[nodiscard]] __host__ __device__ inline auto sinh_elem(M&& m) noexcept requires(CUDA<M>) {
        return device_obj<MatrixExpr<ExprType::Sinh, M&&>>(std::forward<M>(m));
    }
}
