/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.h"

namespace Physica {
    template<Matrix M>
    class MatrixExpr<ExprID::Sin, M>
            : public UnaryMatrixExpr<ExprID::Sin, M> {
        using Base = UnaryMatrixExpr<ExprID::Sin, M>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return sin(Base::getExpr().calc(row, col)); }

        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Matrix M>
    auto MatrixExpr<ExprID::Sin, M>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return sin_elem(std::forward<Self>(self).getExpr().values());
    }

    template<Matrix M>
    [[nodiscard, gnu::always_inline]] auto sin_elem(M&& m) noexcept requires(!DeviceObj<M>) {
        return MatrixExpr<ExprID::Sin, M&&>(std::forward<M>(m));
    }
}
