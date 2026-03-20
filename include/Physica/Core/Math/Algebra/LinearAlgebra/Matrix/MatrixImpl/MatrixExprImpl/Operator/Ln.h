/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.h"

namespace Physica {
    template<Matrix M>
    class MatrixExpr<ExprID::Ln, M>
            : public UnitaryMatrixExpr<ExprID::Ln, M> {
        using Base = UnitaryMatrixExpr<ExprID::Ln, M>;
    public:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return ln(Base::getExpr().calc(row, col)); }

        [[nodiscard]] Tv calc_value(size_t row, size_t col) const {
            return ln(Base::getExpr().calc_value(row, col));
        }

        void reverse(const Vector auto& grad) const noexcept;
    };

    template<Matrix M>
    void MatrixExpr<ExprID::Ln, M>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& expr = Base::getExpr();
        expr.reverse(divide(grad, expr.values()));
    }

    template<Matrix M>
    [[nodiscard]] auto ln_elem(M&& m) noexcept requires(!DeviceObj<M>) {
        return MatrixExpr<ExprID::Ln, M&&>(std::forward<M>(m));
    }
}
