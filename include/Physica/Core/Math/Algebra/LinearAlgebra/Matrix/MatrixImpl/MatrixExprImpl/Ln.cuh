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

#include "../MatrixExpr.h"

namespace Physica {
    template<Matrix M>
    class device_obj<MatrixExpr<ExprType::Ln, M>>
            : public device_obj<UnitaryMatrixExpr<ExprType::Ln, M>> {
        using Base = device_obj<UnitaryMatrixExpr<ExprType::Ln, M>>;
    public:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const { return ln(Base::getExpr().calc(row, col)); }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            return ln(Base::getExpr().calc_value(row, col));
        }

        template<Vector U>
        void reverse(const U& grad) const noexcept requires(isReverseDiff);
    };

    template<Matrix M>
    template<Vector U>
    void device_obj<MatrixExpr<ExprType::Ln, M>>::reverse(const U& grad) const noexcept requires(isReverseDiff) {
        const auto& expr = Base::getExpr();
        expr.reverse(divide(grad, expr.values()));
    }

    template<Matrix M>
    [[nodiscard]] inline auto ln_elem(M&& m) noexcept requires(CUDA<M>) {
        return device_obj<MatrixExpr<ExprType::Ln, M&&>>(std::forward<M>(m));
    }
}
