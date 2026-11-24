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
    class MatrixExpr<ExprID::Relu, M>
            : public UnaryMatrixExpr<ExprID::Relu, M> {
        using Base = UnaryMatrixExpr<ExprID::Relu, M>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return relu(Base::getExpr().calc(row, col)); }

        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Matrix M>
    auto MatrixExpr<ExprID::Relu, M>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return relu_elem(std::forward<Self>(self).getExpr().values());
    }


    template<Matrix M>
    [[nodiscard, gnu::always_inline]] auto relu_elem(M&& m) noexcept requires(!DeviceObj<M>) {
        return MatrixExpr<ExprID::Relu, M&&>(std::forward<M>(m));
    }
}
