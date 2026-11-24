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
    class VectorExpr<ExprID::Ln1p, V> : public UnaryVectorExpr<ExprID::Ln1p, V> {
        using This = VectorExpr<ExprID::Ln1p, V>;
        using Base = UnaryVectorExpr<ExprID::Ln1p, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Trv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto input) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input, size_t count) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const { return ln1p(Base::getExpr().calc(index)); }

        void reverse(const Scalar auto& grad) const noexcept;
        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Ln1p, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return ln1p(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Ln1p, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return ln1p(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Ln1p, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return ln1p(input.template load<Size>(count));
    }

    template<Vector V>
    void VectorExpr<ExprID::Ln1p, V>::reverse(const Scalar auto& grad) const noexcept {
        static_assert(Base::isReverseDiff());
        const auto& expr = Base::getExpr();
        expr.reverse(divide(grad.value(), (Trv(1) + expr.values())));
    }

    template<Vector V>
    auto VectorExpr<ExprID::Ln1p, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return ln1p(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto ln1p(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Ln1p, V&&>(std::forward<V>(v));
    }
}
