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
    class VectorExpr<ExprID::Unit, V>
            : public UnaryVectorExpr<ExprID::Unit, V> {
        using Base = UnaryVectorExpr<ExprID::Unit, V>;
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
        [[nodiscard]] T calc(size_t i) const { return unit(Base::getExpr().calc(i)); }
    };

    template<Vector V>
    auto VectorExpr<ExprID::Unit, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return unit(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Unit, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return unit(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Unit, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return unit(input.template load<Size>(count)).cutoff(count);
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto unit(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Unit, V&&>(std::forward<V>(v));
    }
}
