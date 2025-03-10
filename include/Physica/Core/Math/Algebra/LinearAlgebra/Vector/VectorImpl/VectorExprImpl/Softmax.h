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
    class VectorExpr<ExprType::Softmax, V>
            : public UnitaryVectorExpr<ExprType::Softmax, V> {
        using Base = UnitaryVectorExpr<ExprType::Softmax, V>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        template<Vector V1, class Executor = SeqExecutor>
        inline void assign(V1& v) const;

        [[nodiscard]] T calc(size_t i) const { return Base::getExpr().softmax(i); }
        [[nodiscard]] Tv calc_value(size_t i) const { return Base::getExpr().values().softmax(i); }
        [[nodiscard]] T calc(size_t i, T lnsumexp) const;
        [[nodiscard]] Tv calc_value(size_t i, Tv lnsumexp) const;

        template<Vector U, Vector V1>
        void reverse(const U& y, const V1& grad_) const noexcept requires(isReverseDiff) ;
    };

    template<Vector V>
    auto VectorExpr<ExprType::Softmax, V>::calc(size_t i, T lnsumexp) const -> T {
        return exp(Base::getExpr().calc(i) - lnsumexp);
    }

    template<Vector V>
    auto VectorExpr<ExprType::Softmax, V>::calc_value(size_t i, Tv lnsumexp) const -> Tv {
        return exp(Base::getExpr().calc_value(i) - lnsumexp);
    }

    template<Vector V>
    template<Vector V1, class Executor>
    inline void VectorExpr<ExprType::Softmax, V>::assign(V1& v) const {
        const auto& expr = Base::getExpr();
        const T factor = expr.lnSumExp();
        v = exp(expr - factor);
    }

    template<Vector V>
    template<Vector U, Vector V1>
    void VectorExpr<ExprType::Softmax, V>::reverse(const U& y, const V1& grad_) const noexcept requires(isReverseDiff) {
        const auto& grad = grad_.values();
        Base::getExpr().reverse(hadamard(y, grad) - (y * grad) * y);
    }

    template<Vector V>
    [[nodiscard]] inline auto softmax(V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Softmax, V&&>(std::forward<V>(v));
    }
}
