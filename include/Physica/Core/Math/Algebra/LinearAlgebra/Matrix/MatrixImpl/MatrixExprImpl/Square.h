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

#include "../MatrixExpr.h"

namespace Physica {
    template<Matrix M>
    class MatrixExpr<ExprType::Square, M>
            : public UnitaryMatrixExpr<ExprType::Square, M> {
        using Base = UnitaryMatrixExpr<ExprType::Square, M>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept requires(isReverseDiff);
        using Base::reverse;
    };

    template<Matrix M>
    auto MatrixExpr<ExprType::Square, M>::calc(size_t row, size_t col) const -> T {
        return square(Base::getExpr().calc(row, col));
    }

    template<Matrix M>
    auto MatrixExpr<ExprType::Square, M>::calc_value(size_t row, size_t col) const -> Tv {
        return square(Base::getExpr().calc_value(row, col));
    }

    template<Matrix M>
    void MatrixExpr<ExprType::Square, M>::reverse(const Matrix auto& grad) const noexcept requires(isReverseDiff) {
        const auto& expr = Base::getExpr();
        expr.reverse(Tv(2) * hadamard(expr.values(), grad));
    }

    template<Matrix M>
    [[nodiscard]] auto square_elem(M&& m) noexcept requires(!CUDA<M>) {
        return MatrixExpr<ExprType::Square, M&&>(std::forward<M>(m));
    }
}
