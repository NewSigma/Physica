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
    class VectorExpr<ExprID::Ln, V> : public UnaryVectorExpr<ExprID::Ln, V> {
        using This = VectorExpr<ExprID::Ln, V>;
        using Base = UnaryVectorExpr<ExprID::Ln, V>;
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
        [[nodiscard]] CoDiff<T> calc(size_t index) const;

        void reverse(const auto& grad) const noexcept;
        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Ln, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return ln(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Ln, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        auto x = input.template load<Size>();
        assert(x.isPositive().horizontal_and());
        return ln(x);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Ln, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return ln(input.template load<Size>(count)).cutoff(count);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Ln, V>::calc(size_t index) const -> CoDiff<T> {
        return ln(Base::getExpr().calc(index));
    }

    template<Vector V>
    void VectorExpr<ExprID::Ln, V>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& expr = Base::getExpr();
        if constexpr (Scalar<decltype(grad)>)
            expr.reverse(divide(grad.value(), expr.values()));
        else {
            static_assert(Vector<decltype(grad)>, "[Error]: Unexpected type");
            expr.reverse(divide(grad.values(), expr.values()));
        }
    }

    template<Vector V>
    auto VectorExpr<ExprID::Ln, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return ln(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto ln(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Ln, V&&>(std::forward<V>(v));
    }
}
