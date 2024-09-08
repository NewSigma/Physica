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
    template<class MatrixType>
    class MatrixExpr<ExpressionType::Ln, MatrixType>
            : public UnitaryMatrixExpr<ExpressionType::Ln, MatrixType> {
        using Base = UnitaryMatrixExpr<ExpressionType::Ln, MatrixType>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return ln(Base::getExpr().calc(row, col)); }
    };

    template<class MatrixType>
    [[nodiscard]] inline auto ln_elem(const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExpr<ExpressionType::Ln, MatrixType>(m.getDerived());
    }
}
