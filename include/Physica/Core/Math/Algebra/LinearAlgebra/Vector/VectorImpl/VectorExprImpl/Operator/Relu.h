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
    class VectorExpr<ExprID::Relu, V> : public UnaryVectorExpr<ExprID::Relu, V> {
        using This = VectorExpr<ExprID::Relu, V>;
        using Base = UnaryVectorExpr<ExprID::Relu, V>;
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

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        using Base::getExpr;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Relu, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return relu(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Relu, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return relu(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Relu, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return relu(input.template load<Size>(count));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Relu, V>::calc(size_t index) const -> CoDiff<T> {
        return relu(getExpr().calc(index));
    }

    template<Vector V>
    void VectorExpr<ExprID::Relu, V>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& g = grad.values();
        for (size_t i = 0; i < Base::getLength(); ++i) {
            auto v = getExpr().calc(i);
            v.reverse(v.isPositive() ? g.calc(i) : Tv(0));
        }
    }

    template<Vector V>
    auto VectorExpr<ExprID::Relu, V>::values(this auto&& self) noexcept {
        return relu(std::forward<decltype(self)>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto relu(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Relu, V&&>(std::forward<V>(v));
    }
}
