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
    class VectorExpr<ExprID::Cbrt, V> : public UnaryVectorExpr<ExprID::Cbrt, V> {
        using This = VectorExpr<ExprID::Cbrt, V>;
        using Base = UnaryVectorExpr<ExprID::Cbrt, V>;
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
        [[nodiscard]] CoDiff<T> calc(size_t s) const { return cbrt(Base::getExpr().calc(s)); }
        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Cbrt, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return cbrt(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Cbrt, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return cbrt(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Cbrt, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return cbrt(input.template load<Size>(count));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Cbrt, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return cbrt(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto cbrt(V&& v) noexcept {
        return VectorExpr<ExprID::Cbrt, V&&>(std::forward<V>(v));
    }
}
