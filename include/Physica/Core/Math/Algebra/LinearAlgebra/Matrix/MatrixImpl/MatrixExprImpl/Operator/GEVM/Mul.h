/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../Mul.h"

namespace Physica {
    template<Matrix M, Scalar U> requires(instanceof<M, GEVM>)
    class MatrixExpr<ExprID::Mul, M, U> : public BinaryMatrixExpr<ExprID::Mul, M, U> {
        using Base = BinaryMatrixExpr<ExprID::Mul, M, U>;
        using This = MatrixExpr<ExprID::Mul, M, U>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_mkl(Matrix auto& target) const;
        void assign_base(Matrix auto& target) const;
        void assign_add(Matrix auto& target) const;

        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;

        [[nodiscard]] auto values(this auto&& self) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M, Scalar U> requires(instanceof<M, GEVM>)
    void MatrixExpr<ExprID::Mul, M, U>::assign(Matrix auto& target) const {
        assign_base(target);
    }

    template<Matrix M, Scalar U> requires(instanceof<M, GEVM>)
    void MatrixExpr<ExprID::Mul, M, U>::assign_base(Matrix auto& target) const {
        const auto& vec = getLHS().getLHS();
        const auto& mat = getLHS().getRHS();
        ((getRHS() * vec) * mat).assign(target);
    }

    template<Matrix M, Scalar U> requires(instanceof<M, GEVM>)
    void MatrixExpr<ExprID::Mul, M, U>::assign_add(Matrix auto& target) const {
        const auto& vec = getLHS().getLHS();
        const auto& mat = getLHS().getRHS();
        ((getRHS() * vec) * mat).assign_add(target);
    }

    template<Matrix M, Scalar U> requires(instanceof<M, GEVM>)
    auto MatrixExpr<ExprID::Mul, M, U>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return getLHS().calc(row, col) * getRHS();
    }

    template<Matrix M, Scalar U> requires(instanceof<M, GEVM>)
    auto MatrixExpr<ExprID::Mul, M, U>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().value();
    }
}
