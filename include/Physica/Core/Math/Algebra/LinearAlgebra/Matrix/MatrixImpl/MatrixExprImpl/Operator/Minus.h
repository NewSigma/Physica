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
    class MatrixExpr<ExprID::Minus, M>
            : public UnaryMatrixExpr<ExprID::Minus, M> {
        using Base = UnaryMatrixExpr<ExprID::Minus, M>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return -Base::getExpr().calc(row, col); }

        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Matrix M>
    auto MatrixExpr<ExprID::Minus, M>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return -std::forward<Self>(self).getExpr().values();
    }
}
