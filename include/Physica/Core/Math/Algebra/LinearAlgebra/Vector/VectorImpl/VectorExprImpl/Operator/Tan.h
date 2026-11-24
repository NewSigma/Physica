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
    class VectorExpr<ExprID::Tan, V> : public UnaryVectorExpr<ExprID::Tan, V> {
        using This = VectorExpr<ExprID::Tan, V>;
        using Base = UnaryVectorExpr<ExprID::Tan, V>;
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
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return tan(Base::getExpr().calc(index)); }
        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Tan, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return tan(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Tan, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return tan(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Tan, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return tan(input.template load<Size>(count));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Tan, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return tan(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto tan(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Tan, V&&>(std::forward<V>(v));
    }
}
