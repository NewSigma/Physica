/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.h"

namespace Physica {
    template<Vector V>
    class VectorExpr<ExprID::Softmax, V>
            : public UnaryVectorExpr<ExprID::Softmax, V> {
        using Base = UnaryVectorExpr<ExprID::Softmax, V>;
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

        [[nodiscard]] CoDiff<T> calc(size_t i) const { return Base::getExpr().softmax(i); }
        [[nodiscard]] Tv calc_value(size_t i) const { return Base::calc_value(i); }
        [[nodiscard]] CoDiff<T> calc(size_t i, T lnsumexp) const;
        [[nodiscard]] Tv calc_value(size_t i, Tv lnsumexp) const;

        using Base::reverse;
        void reverse(const Vector auto& y, const Vector auto& grad) const noexcept;
        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Static members */
        __host__ __device__ consteval static bool isFastPacket() noexcept { return false; }
    };

    template<Vector V>
    auto VectorExpr<ExprID::Softmax, V>::calc(size_t i, T lnsumexp) const -> CoDiff<T> {
        return exp(Base::getExpr().calc(i) - lnsumexp);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Softmax, V>::calc_value(size_t i, Tv lnsumexp) const -> Tv {
        return values().calc(i, lnsumexp);
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
        static_assert(isReverseDiff());
        const auto& g = grad.values();
        Base::getExpr().reverse(hadamard(y, g) - (y * g) * y);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Softmax, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return softmax(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto softmax(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Softmax, V&&>(std::forward<V>(v));
    }
}
