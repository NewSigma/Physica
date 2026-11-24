/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.h"

namespace Physica {
    template<Matrix M>
    class MatrixExpr<ExprID::Unit, M>
            : public UnaryMatrixExpr<ExprID::Unit, M> {
        using Base = UnaryMatrixExpr<ExprID::Unit, M>;
    public:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const {
            return unit(Base::getExpr().calc(row, col));
        }

        [[nodiscard]] auto values(this auto&& self) noexcept;
    };

    template<Matrix M>
    auto MatrixExpr<ExprID::Unit, M>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return unit_elem(std::forward<Self>(self).getExpr().values());
    }

    template<Matrix M>
    [[nodiscard, gnu::always_inline]] auto unit_elem(M&& m) noexcept requires(!DeviceObj<M>) {
        return MatrixExpr<ExprID::Unit, M&&>(std::forward<M>(m));
    }
}
