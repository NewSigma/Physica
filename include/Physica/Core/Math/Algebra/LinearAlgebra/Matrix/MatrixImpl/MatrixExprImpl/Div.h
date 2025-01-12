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

namespace Physica::Core {
    template<class T, class U>
    class MatrixExpr<ExprType::Div, T, U>
            : public BinaryMatrixExpr<ExprType::Div, T, U> {
        static_assert(Scalar<T> || Scalar<U>, "[Error]: Either types should be Scalar");

        using Base = BinaryMatrixExpr<ExprType::Div, T, U>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    public:
        MatrixExpr(const T& lhs, const U& rhs) : Base(lhs, rhs) {
            if constexpr (Matrix<T>)
                assert(!rhs.isZero() && "[Error]: Divide by zero");
        }
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc(row, col) / Base::getRHS();
            else
                return Base::getLHS() / Base::getRHS().calc(row, col);
        }

        [[nodiscard]] ValueType calc_value(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc_value(row, col) / Base::getRHS().value();
            else
                return Base::getLHS().value() / Base::getRHS().calc_value(row, col);
        }
    };

    template<Matrix T, Scalar U>
    [[nodiscard]] inline auto operator/(const T& m, const U& x) noexcept {
        return MatrixExpr<ExprType::Div, T, U>(m, x);
    }

    template<Scalar T, Matrix U>
    [[nodiscard]] inline auto operator/(const T& m, const U& x) noexcept {
        return MatrixExpr<ExprType::Div, T, U>(m, x);
    }
}
