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
    class VectorExpr<ExprID::LnCosh, V> : public UnaryVectorExpr<ExprID::LnCosh, V> {
        using This = VectorExpr<ExprID::LnCosh, V>;
        using Base = UnaryVectorExpr<ExprID::LnCosh, V>;
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

        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::LnCosh, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return lncosh(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::LnCosh, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return lncosh(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::LnCosh, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return lncosh(input.template load<Size>(count));
    }

    template<Vector V>
    auto VectorExpr<ExprID::LnCosh, V>::calc(size_t index) const -> CoDiff<T> {
        return lncosh(Base::getExpr().calc(index));
    }

    template<Vector V>
    auto VectorExpr<ExprID::LnCosh, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return lncosh(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto lncosh(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::LnCosh, V&&>(std::forward<V>(v));
    }
}
