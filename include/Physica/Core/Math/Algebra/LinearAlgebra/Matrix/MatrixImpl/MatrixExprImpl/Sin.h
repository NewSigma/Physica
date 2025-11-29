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
    class MatrixExpr<ExprID::Sin, M>
            : public UnitaryMatrixExpr<ExprID::Sin, M> {
        using Base = UnitaryMatrixExpr<ExprID::Sin, M>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return sin(Base::getExpr().calc(row, col)); }

        [[nodiscard]] Tv calc_value(size_t row, size_t col) const {
            return sin(Base::getExpr().calc_value(row, col));
        }
    };

    template<Matrix M>
    [[nodiscard]] auto sin_elem(M&& m) noexcept requires(!CUDA<M>) {
        return MatrixExpr<ExprID::Sin, M&&>(std::forward<M>(m));
    }
}
