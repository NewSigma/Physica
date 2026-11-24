/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.h"

namespace Physica {
    template<Matrix M>
    class MatrixExpr<ExprID::Square, M>
            : public UnaryMatrixExpr<ExprID::Square, M> {
        using Base = UnaryMatrixExpr<ExprID::Square, M>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Matrix M>
    auto MatrixExpr<ExprID::Square, M>::calc(size_t row, size_t col) const -> T {
        return square(Base::getExpr().calc(row, col));
    }

    template<Matrix M>
    void MatrixExpr<ExprID::Square, M>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& expr = Base::getExpr();
        expr.reverse(Tv(2) * hadamard(expr.values(), grad));
    }

    template<Matrix M>
    auto MatrixExpr<ExprID::Square, M>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return square_elem(std::forward<Self>(self).getExpr().values());
    }


    template<Matrix M>
    [[nodiscard, gnu::always_inline]] auto square_elem(M&& m) noexcept requires(!DeviceObj<M>) {
        return MatrixExpr<ExprID::Square, M&&>(std::forward<M>(m));
    }
}
