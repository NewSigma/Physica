/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.h"

namespace Physica {
    template<Vector V>
    class VectorExpr<ExprID::ArcSinh, V> : public UnaryVectorExpr<ExprID::ArcSinh, V> {
        using This = VectorExpr<ExprID::ArcSinh, V>;
        using Base = UnaryVectorExpr<ExprID::ArcSinh, V>;
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
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return arcsinh(Base::getExpr().calc(index)); }

        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::ArcSinh, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return arcsinh(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::ArcSinh, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return arcsinh(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::ArcSinh, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return arcsinh(input.template load<Size>(count));
    }

    template<Vector V>
    auto VectorExpr<ExprID::ArcSinh, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return arcsinh(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto arcsinh(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::ArcSinh, V&&>(std::forward<V>(v));
    }
}
