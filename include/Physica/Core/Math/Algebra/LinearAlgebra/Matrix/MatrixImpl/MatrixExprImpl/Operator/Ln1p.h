/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
        using Base::isReverseDiff;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;

        void reverse(const auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] auto row(this auto&&, size_t r) noexcept;
        [[nodiscard]] auto col(this auto&&, size_t c) noexcept;

        [[nodiscard]] CoDiff<T> sum() const;

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
    auto MatrixExpr<ExprID::Ln1p, M>::sum() const -> CoDiff<T> {
        if constexpr (isReverseDiff())
            return Base::sum();
        else {
            T result = 0;
            if constexpr (This::isColMatrix()) {
                for (size_t c = 0; c < Base::getCol(); ++c)
                    result += col(c).sum();
            }
            else {
                for (size_t r = 0; r < Base::getRow(); ++r)
                    result += row(r).sum();
            }
            return result;
        }
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
