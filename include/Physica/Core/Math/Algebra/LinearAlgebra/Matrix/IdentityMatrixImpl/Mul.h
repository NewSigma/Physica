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

#include "../IdentityMatrix.h"

namespace Physica {
    template<Matrix M, Scalar U> requires(instanceof_tx<IdentityMatrix, M>)
    class MatrixExpr<ExprID::Mul, M, U>
            : public BinaryMatrixExpr<ExprID::Mul, M, U> {
        using Base = BinaryMatrixExpr<ExprID::Mul, M, U>;
        using This = MatrixExpr<ExprID::Mul, M, U>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        using Base::operator*;
        [[nodiscard]] auto operator*(Scalar auto x) const noexcept;
        [[nodiscard]] auto operator-() const& noexcept;
        [[nodiscard]] auto operator-() && noexcept;
        /* Operations */
        void assign(Matrix auto&& target) const;
        void assign_add(Matrix auto&& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M, Scalar U> requires(instanceof_tx<IdentityMatrix, M>)
    auto MatrixExpr<ExprID::Mul, M, U>::operator*(Scalar auto x) const noexcept {
        return getLHS() * (x * getRHS());
    }

    template<Matrix M, Scalar U> requires(instanceof_tx<IdentityMatrix, M>)
    auto MatrixExpr<ExprID::Mul, M, U>::operator-() const& noexcept {
        return getLHS() * (-getRHS());
    }

    template<Matrix M, Scalar U> requires(instanceof_tx<IdentityMatrix, M>)
    auto MatrixExpr<ExprID::Mul, M, U>::operator-() && noexcept {
        return std::move(getLHS()) * (-getRHS());
    }

    template<Matrix M, Scalar U> requires(instanceof_tx<IdentityMatrix, M>)
    void MatrixExpr<ExprID::Mul, M, U>::assign(Matrix auto&& target) const {
        target.assert_assign(*this);
        target.zeros();
        assign_add(target);
    }

    template<Matrix M, Scalar U> requires(instanceof_tx<IdentityMatrix, M>)
    void MatrixExpr<ExprID::Mul, M, U>::assign_add(Matrix auto&& target) const {
        target.assert_assign(*this);
        target.diag() += getRHS();
    }

    template<Matrix M, Scalar U> requires(instanceof_tx<IdentityMatrix, M>)
    auto MatrixExpr<ExprID::Mul, M, U>::calc(size_t row, size_t col) const -> T {
        return row != col ? T(0) : getRHS();
    }

    template<Matrix M, Scalar U> requires(instanceof_tx<IdentityMatrix, M>)
    auto MatrixExpr<ExprID::Mul, M, U>::calc_value(size_t row, size_t col) const -> Tv {
        return row != col ? Tv(0) : getRHS().value();
    }
}
