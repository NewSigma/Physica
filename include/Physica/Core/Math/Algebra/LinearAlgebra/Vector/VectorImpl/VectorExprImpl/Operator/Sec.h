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
    class VectorExpr<ExprID::Sec, V> : public UnaryVectorExpr<ExprID::Sec, V> {
        using This = VectorExpr<ExprID::Sec, V>;
        using Base = UnaryVectorExpr<ExprID::Sec, V>;
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
        [[nodiscard]] CoDiff<T> calc(size_t index) const {
            return sec(Base::getExpr().calc(index));
        }

        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Sec, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return sec(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Sec, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return sec(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Sec, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return sec(input.template load<Size>(count)).cutoff(count);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Sec, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return sec(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto sec(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Sec, V&&>(std::forward<V>(v));
    }
}
