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
    class VectorExpr<ExprID::Cos, V> : public UnaryVectorExpr<ExprID::Cos, V> {
        using This = VectorExpr<ExprID::Cos, V>;
        using Base = UnaryVectorExpr<ExprID::Cos, V>;
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
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return cos(Base::getExpr().calc(index)); }
        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Cos, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return cos(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Cos, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return cos(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Cos, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return cos(input.template load<Size>(count)).cutoff(count);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Cos, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return cos(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto cos(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Cos, V&&>(std::forward<V>(v));
    }
}
