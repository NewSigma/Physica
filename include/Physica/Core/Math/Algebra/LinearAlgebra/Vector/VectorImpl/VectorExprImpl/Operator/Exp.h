/*
 * Copyright 2024-2026 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.h"

namespace Physica {
    template<Vector V>
    class VectorExpr<ExprID::Exp, V> : public UnitaryVectorExpr<ExprID::Exp, V> {
        using This = VectorExpr<ExprID::Exp, V>;
        using Base = UnitaryVectorExpr<ExprID::Exp, V>;
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
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const noexcept;
        void assign_mkl(Vector auto& v) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;

        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Exp, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return exp(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Exp, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        return exp(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Exp, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        return exp(input.template load<Size>(count)).cutoff(count);
    }

    template<Vector V>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Exp, V>::assign(Vector auto&& v) const noexcept {
        using V1 = std::remove_cvref_t<decltype(v)>;
        constexpr size_t Length = std::max(Base::SizeAtCompile, V1::SizeAtCompile);
        constexpr bool SmallVector = 0 < Length && Length <= 32;
        if constexpr (Internal::EnableMKL<V, V1>::value && !SmallVector && T::Prec == Float64) {
            if (Base::getLength() <= 32)
                Base::template assign_base<P>(v);
            else
                assign_mkl(v);
        }
        else
            Base::template assign_base<P>(v);
    }

    template<Vector V>
    auto VectorExpr<ExprID::Exp, V>::calc(size_t index) const -> CoDiff<T> {
        return exp(Base::getExpr().calc(index));
    }

    template<Vector V>
    void VectorExpr<ExprID::Exp, V>::reverse(const auto& grad) const noexcept {
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
    auto VectorExpr<ExprID::Exp, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return exp(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto exp(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Exp, V&&>(std::forward<V>(v));
    }
}

#ifdef PHYSICA_MKL
    #include "MKL/Exp.h"
#endif
