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
    class VectorExpr<ExprID::Sin, V> : public UnaryVectorExpr<ExprID::Sin, V> {
        using This = VectorExpr<ExprID::Sin, V>;
        using Base = UnaryVectorExpr<ExprID::Sin, V>;
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
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return sin(Base::getExpr().calc(index)); }
        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Sin, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return sin(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Sin, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return sin(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Sin, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return sin(input.template load<Size>(count));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Sin, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return sin(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto sin(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Sin, V&&>(std::forward<V>(v));
    }
}
