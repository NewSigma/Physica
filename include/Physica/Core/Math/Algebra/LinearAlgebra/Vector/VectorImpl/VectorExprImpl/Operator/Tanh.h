/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.h"

namespace Physica {
    template<Vector V>
    class VectorExpr<ExprID::Tanh, V> : public UnaryVectorExpr<ExprID::Tanh, V> {
        using This = VectorExpr<ExprID::Tanh, V>;
        using Base = UnaryVectorExpr<ExprID::Tanh, V>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto input) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input, size_t count) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return tanh(Base::getExpr().calc(index)); }

        using Base::reverse;
        void reverse(const Vector auto& y, const Vector auto& grad) const noexcept;
        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Tanh, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return tanh(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Tanh, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return tanh(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Tanh, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return tanh(input.template load<Size>(count));
    }

    template<Vector V>
    void VectorExpr<ExprID::Tanh, V>::reverse(const Vector auto& y, const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& expr = Base::getExpr();
        expr.reverse(hadamard((Tv(1) - square(y.values())), grad));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Tanh, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return tanh(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto tanh(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Tanh, V&&>(std::forward<V>(v));
    }
}
