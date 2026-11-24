/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/TensorImpl/TensorExpr.h"

namespace Physica {
    template<Tensor X1, Tensor X2>
    class TensorExpr<ExprID::Add, X1, X2>
            : public BinaryTensorExpr<ExprID::Add, X1, X2> {
        static_assert(canonicalized<X1, X2>());
        using Base = BinaryTensorExpr<ExprID::Add, X1, X2>;
    public:
        using typename Base::IndexType;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] T calc(const IndexType& indices) const {
            return Base::getLHS().calc(indices) + Base::getRHS().calc(indices);
        }
    };

    template<Tensor X, Tensor Y>
    [[nodiscard, gnu::always_inline]] auto operator+(X&& x, Y&& y) noexcept requires(!DeviceObj<X> && !DeviceObj<Y>) {
        if constexpr (!canonicalized(x, y))
            return std::forward<Y>(y) + std::forward<X>(x);
        else
            return TensorExpr<ExprID::Add, X&&, Y&&>(std::forward<X>(x), std::forward<Y>(y));
    }
}
