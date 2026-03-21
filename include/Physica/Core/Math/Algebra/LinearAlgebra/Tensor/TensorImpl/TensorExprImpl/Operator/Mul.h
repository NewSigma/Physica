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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/TensorImpl/TensorExpr.h"

namespace Physica {
    template<Tensor X, Scalar U>
    class TensorExpr<ExprID::Mul, X, U>
            : public BinaryTensorExpr<ExprID::Mul, X, U> {
        using Base = BinaryTensorExpr<ExprID::Mul, X, U>;
    public:
        using typename Base::IndexType;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] T calc(const IndexType& indices) const {
            return Base::getLHS().calc(indices) * Base::getRHS();
        }
    };

    template<Tensor X1, Tensor X2>
    class TensorExpr<ExprID::Mul, X1, X2>
            : public BinaryTensorExpr<ExprID::Mul, X1, X2> {
        using Base = BinaryTensorExpr<ExprID::Mul, X1, X2>;
    public:
        using typename Base::IndexType;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] T calc(const IndexType& indices) const {
            return Base::getLHS().calc(indices) * Base::getRHS().calc(indices);
        }
    };

    template<Tensor X, Scalar U>
    [[nodiscard]] auto operator*(X&& x, U&& y) noexcept {
        return TensorExpr<ExprID::Mul, X&&, U&&>(std::forward<X>(x), std::forward<U>(y));
    }

    template<Tensor X, Scalar U>
    [[nodiscard]] auto operator*(U&& y, X&& x) noexcept {
        return std::forward<X>(x) * std::forward<U>(y);
    }

    template<Tensor X, Tensor Y>
    [[nodiscard]] auto hadamard(X&& x, Y&& y) noexcept {
        if constexpr (!canonicalized(x, y))
            return hadamard(std::forward<Y>(y), std::forward<X>(x));
        else
            return TensorExpr<ExprID::Mul, X&&, Y&&>(std::forward<X>(x), std::forward<Y>(y));
    }
}
