/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.h"

namespace Physica {
    template<Vector V>
    class VectorExpr<ExprID::Sigmoid, V> : public UnaryVectorExpr<ExprID::Sigmoid, V> {
        using This = VectorExpr<ExprID::Sigmoid, V>;
        using Base = UnaryVectorExpr<ExprID::Sigmoid, V>;
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
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return sigmoid(Base::getExpr().calc(index)); }
        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Sigmoid, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return sigmoid(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Sigmoid, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return sigmoid(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Sigmoid, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return sigmoid(input.template load<Size>(count)).cutoff(count);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Sigmoid, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return sigmoid(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto sigmoid(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Sigmoid, V&&>(std::forward<V>(v));
    }
}
