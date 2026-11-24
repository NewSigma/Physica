/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../IdentityMatrix.h"

namespace Physica {
    template<Matrix M, Scalar U> requires(instanceof_tx<M, IdentityMatrix>)
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
        [[nodiscard]] auto operator-(this auto&&) noexcept;
        /* Operations */
        void assign(Matrix auto&& target) const;
        void assign_add(Matrix auto&& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;

        [[nodiscard]] auto values(this auto&& self) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M, Scalar U> requires(instanceof_tx<M, IdentityMatrix>)
    auto MatrixExpr<ExprID::Mul, M, U>::operator*(Scalar auto x) const noexcept {
        return getLHS() * (x * getRHS());
    }

    template<Matrix M, Scalar U> requires(instanceof_tx<M, IdentityMatrix>)
    auto MatrixExpr<ExprID::Mul, M, U>::operator-(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS() * (-std::forward<Self>(self).getRHS());
    }

    template<Matrix M, Scalar U> requires(instanceof_tx<M, IdentityMatrix>)
    void MatrixExpr<ExprID::Mul, M, U>::assign(Matrix auto&& target) const {
        target.assert_assign(*this);
        target.zeros();
        assign_add(target);
    }

    template<Matrix M, Scalar U> requires(instanceof_tx<M, IdentityMatrix>)
    void MatrixExpr<ExprID::Mul, M, U>::assign_add(Matrix auto&& target) const {
        target.assert_assign(*this);
        target.diag() += getRHS();
    }

    template<Matrix M, Scalar U> requires(instanceof_tx<M, IdentityMatrix>)
    auto MatrixExpr<ExprID::Mul, M, U>::calc(size_t row, size_t col) const -> T {
        return row != col ? T(0) : getRHS();
    }

    template<Matrix M, Scalar U> requires(instanceof_tx<M, IdentityMatrix>)
    auto MatrixExpr<ExprID::Mul, M, U>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }
}
