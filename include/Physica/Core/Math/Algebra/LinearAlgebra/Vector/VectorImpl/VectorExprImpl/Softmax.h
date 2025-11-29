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

#include "../VectorExpr.h"

namespace Physica {
    template<Vector V>
    class VectorExpr<ExprID::Softmax, V>
            : public UnitaryVectorExpr<ExprID::Softmax, V> {
        using Base = UnitaryVectorExpr<ExprID::Softmax, V>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const;

        [[nodiscard]] T calc(size_t i) const { return Base::getExpr().softmax(i); }
        [[nodiscard]] Tv calc_value(size_t i) const { return Base::getExpr().values().softmax(i); }
        [[nodiscard]] T calc(size_t i, T lnsumexp) const;
        [[nodiscard]] Tv calc_value(size_t i, Tv lnsumexp) const;

        void reverse(const Vector auto& y, const Vector auto& grad) const noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Softmax, V>::calc(size_t i, T lnsumexp) const -> T {
        return exp(Base::getExpr().calc(i) - lnsumexp);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Softmax, V>::calc_value(size_t i, Tv lnsumexp) const -> Tv {
        return exp(Base::getExpr().calc_value(i) - lnsumexp);
    }

    template<Vector V>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Softmax, V>::assign(Vector auto&& v) const {
        const auto& expr = Base::getExpr();
        const T factor = expr.lnSumExp();
        v = exp(expr - factor);
    }

    template<Vector V>
    void VectorExpr<ExprID::Softmax, V>::reverse(const Vector auto& y, const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& g = grad.values();
        Base::getExpr().reverse(hadamard(y, g) - (y * g) * y);
    }

    template<Vector V>
    [[nodiscard]] auto softmax(V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprID::Softmax, V&&>(std::forward<V>(v));
    }
}
