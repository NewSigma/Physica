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

#include "../MatrixExpr.h"

namespace Physica {
    template<Matrix M>
    class MatrixExpr<ExprType::Ln1p, M>
            : public UnitaryMatrixExpr<ExprType::Ln1p, M> {
        using This = MatrixExpr<ExprType::Ln1p, M>;
        using Base = UnitaryMatrixExpr<ExprType::Ln1p, M>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Trv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        void reverse(const auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] auto row(size_t r) noexcept;
        [[nodiscard]] const auto row(size_t r) const noexcept;
        [[nodiscard]] auto col(size_t c) noexcept;
        [[nodiscard]] const auto col(size_t c) const noexcept;

        [[nodiscard]] auto values() const noexcept;
        /* Getters */
        using Base::getExpr;
    };

    template<Matrix M>
    auto MatrixExpr<ExprType::Ln1p, M>::calc(size_t row, size_t col) const -> T {
        return ln1p(getExpr().calc(row, col));
    }

    template<Matrix M>
    auto MatrixExpr<ExprType::Ln1p, M>::calc_value(size_t row, size_t col) const -> Tv {
        return ln1p(getExpr().calc_value(row, col));
    }

    template<Matrix M>
    void MatrixExpr<ExprType::Ln1p, M>::reverse(const auto& grad) const noexcept {
        static_assert(Base::isReverseDiff);
        using U = decltype(grad);
        auto& expr = Base::getExpr();
        if constexpr (Scalar<U>)
            expr.reverse(reciprocal_elem(expr.values() + Trv(1)) * grad);
        else {
            static_assert(Matrix<U>, "[Error]: Unexpected type");
            expr.reverse(divide_elem(grad, expr.values() + Trv(1)));
        }
    }

    template<Matrix M>
    auto MatrixExpr<ExprType::Ln1p, M>::row(size_t r) noexcept {
        return ln1p(getExpr().row(r));
    }

    template<Matrix M>
    const auto MatrixExpr<ExprType::Ln1p, M>::row(size_t r) const noexcept {
        return const_cast<This&>(*this).row(r);
    }

    template<Matrix M>
    auto MatrixExpr<ExprType::Ln1p, M>::col(size_t c) noexcept {
        return ln1p(getExpr().col(c));
    }

    template<Matrix M>
    const auto MatrixExpr<ExprType::Ln1p, M>::col(size_t c) const noexcept {
        return const_cast<This&>(*this).col(c);
    }

    template<Matrix M>
    auto MatrixExpr<ExprType::Ln1p, M>::values() const noexcept {
        return ln1p_elem(getExpr().values());
    }

    template<Matrix M>
    [[nodiscard]] auto ln1p_elem(M&& m) noexcept requires(!CUDA<M>) {
        return MatrixExpr<ExprType::Ln1p, M&&>(std::forward<M>(m));
    }
}
