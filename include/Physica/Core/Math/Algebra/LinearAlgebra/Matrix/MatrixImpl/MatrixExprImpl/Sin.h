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

namespace Physica {
    template<Matrix T>
    class MatrixExpr<ExprType::Sin, T>
            : public UnitaryMatrixExpr<ExprType::Sin, T> {
        using Base = UnitaryMatrixExpr<ExprType::Sin, T>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return sin(Base::getExpr().calc(row, col)); }

        [[nodiscard]] ValueType calc_value(size_t row, size_t col) const {
            return sin(Base::getExpr().calc_value(row, col));
        }
    };

    template<Matrix T>
    [[nodiscard]] inline auto sin_elem(const T& m) noexcept {
        return MatrixExpr<ExprType::Sin, T>(m);
    }
}
