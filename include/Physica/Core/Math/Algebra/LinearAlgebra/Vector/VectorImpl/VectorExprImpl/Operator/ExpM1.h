/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.h"

namespace Physica {
    template<Vector V>
    class VectorExpr<ExprID::ExpM1, V> : public UnaryVectorExpr<ExprID::ExpM1, V> {
        using This = VectorExpr<ExprID::ExpM1, V>;
        using Base = UnaryVectorExpr<ExprID::ExpM1, V>;
    public:
        using Base::isComplex;
        using Base::isReverseDiff;
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
        void assign_mkl(Vector auto& v) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;

        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::ExpM1, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return expm1(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::ExpM1, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return expm1(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::ExpM1, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return expm1(input.template load<Size>(count));
    }

    template<Vector V>
    auto VectorExpr<ExprID::ExpM1, V>::calc(size_t index) const -> CoDiff<T> {
        return expm1(Base::getExpr().calc(index));
    }

    template<Vector V>
    void VectorExpr<ExprID::ExpM1, V>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff());
        using U = decltype(grad);
        const auto& expr = Base::getExpr();
        if constexpr (Scalar<U>)
            expr.reverse(exp(expr.values()) * grad.value());
        else {
            static_assert(Vector<U>, "[Error]: Unexpected type");
            expr.reverse(hadamard(exp(expr.values()), grad.values()));
        }
    }

    template<Vector V>
    auto VectorExpr<ExprID::ExpM1, V>::values(this auto&& self) noexcept {
        return expm1(std::forward<decltype(self)>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto expm1(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::ExpM1, V&&>(std::forward<V>(v));
    }
}

#ifdef PHYSICA_MKL
    #include "MKL/ExpM1.h"
#endif
