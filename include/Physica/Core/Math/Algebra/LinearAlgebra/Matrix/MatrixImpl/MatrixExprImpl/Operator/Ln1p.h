/*
 * Copyright 2025-2026 Weibo He.
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
    class MatrixExpr<ExprID::Ln1p, M>
            : public UnaryMatrixExpr<ExprID::Ln1p, M> {
        using This = MatrixExpr<ExprID::Ln1p, M>;
        using Base = UnaryMatrixExpr<ExprID::Ln1p, M>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Trv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;

        void reverse(const auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] auto row(this auto&&, size_t r) noexcept;
        [[nodiscard]] auto col(this auto&&, size_t c) noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        using Base::getExpr;
    };

    template<Matrix M>
    auto MatrixExpr<ExprID::Ln1p, M>::calc(size_t row, size_t col) const -> T {
        return ln1p(getExpr().calc(row, col));
    }

    template<Matrix M>
    void MatrixExpr<ExprID::Ln1p, M>::reverse(const auto& grad) const noexcept {
        using U = decltype(grad);
        static_assert(Base::isReverseDiff());
        static_assert(Scalar<U> || Matrix<U>, "[Error]: Unexpected type");
        auto& expr = Base::getExpr();
        expr.reverse(divide_elem(grad, (expr.values() + Trv(1))));
    }

    template<Matrix M>
    auto MatrixExpr<ExprID::Ln1p, M>::row(this auto&& self, size_t r) noexcept {
        using Self = decltype(self);
        return ln1p(std::forward<Self>(self).getExpr().row(r));
    }

    template<Matrix M>
    auto MatrixExpr<ExprID::Ln1p, M>::col(this auto&& self, size_t c) noexcept {
        using Self = decltype(self);
        return ln1p(std::forward<Self>(self).getExpr().col(c));
    }

    template<Matrix M>
    auto MatrixExpr<ExprID::Ln1p, M>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return ln1p_elem(std::forward<Self>(self).getExpr().values());
    }

    template<Matrix M>
    [[nodiscard, gnu::always_inline]] auto ln1p_elem(M&& m) noexcept requires(!DeviceObj<M>) {
        return MatrixExpr<ExprID::Ln1p, M&&>(std::forward<M>(m));
    }
}
