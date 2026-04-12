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
    template<Vector V1, Vector V2>
    class VectorExpr<ExprID::Sub, V1, V2>
            : public BinaryVectorExpr<ExprID::Sub, V1, V2> {
        using Base = BinaryVectorExpr<ExprID::Sub, V1, V2>;
    public:
        using Base::isComplex;
    protected:
        using typename Base::T;
        using typename Base::Tc;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs, size_t count) noexcept;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const;
        void assign_mkl(Vector auto& v) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t s) const;

        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Sub, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept -> CoDiff<T> {
        return *lhs - *rhs;
    }

    template<Vector V1, Vector V2>
    template<int Size>
    auto VectorExpr<ExprID::Sub, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept -> SIMD<T, Size> {
        return lhs.template load<Size>() - rhs.template load<Size>();
    }

    template<Vector V1, Vector V2>
    template<int Size>
    auto VectorExpr<ExprID::Sub, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs, size_t count) noexcept -> SIMD<T, Size> {
        return lhs.template load<Size>(count) - rhs.template load<Size>(count);
    }

    template<Vector V1, Vector V2>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Sub, V1, V2>::assign(Vector auto&& v) const {
        constexpr bool FastAssign1 = Traits<std::remove_cvref_t<V1>>::FastAssign;
        constexpr bool FastAssign2 = Traits<std::remove_cvref_t<V2>>::FastAssign;
        if constexpr (FastAssign1) {
            getLHS().template assign<P>(v);
            v -= getRHS();
        }
        else if constexpr (FastAssign2) {
            static_assert(Traits<decltype(-getRHS())>::FastAssign, "[Debug]: Fast minus implementation is missing");
            (-getRHS()).template assign<P>(v);
            v += getLHS();
        }
        else {
            using V = std::remove_cvref<decltype(v)>::type;
            constexpr size_t SizeAtCompile = Base::getSizeAtCompile(v);
            constexpr size_t Critical = HostDevAttr::LineSizeL1D / sizeof(T);
            constexpr bool UseMKL1 = Internal::EnableMKL<V1, V>::value;
            constexpr bool UseMKL2 = Internal::EnableMKL<V2, V>::value;
            constexpr bool UseMKL3 = SizeAtCompile == Dynamic || SizeAtCompile > Critical;
            constexpr bool UseMKL = UseMKL1 && UseMKL2 && UseMKL3 && (T::Prec != Float64);
            if constexpr (UseMKL) {
                if (Base::getLength() > Critical)
                    assign_mkl(v);
                else
                    Base::template assign_base<P>(v);
            }
            else
                Base::template assign_base<P>(v);
        }
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Sub, V1, V2>::calc(size_t s) const -> CoDiff<T> {
        return getLHS().calc(s) - getRHS().calc(s);
    }

    template<Vector V1, Vector V2>
    void VectorExpr<ExprID::Sub, V1, V2>::reverse(const Vector auto& grad) const noexcept {
        static_assert(Base::isReverseDiff());
        const auto& g = grad.values();
        if constexpr (ReverseDiff<V1>)
            Base::getLHS().reverse(g);
        if constexpr (ReverseDiff<V2>)
            Base::getRHS().reverse(g);
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Sub, V1, V2>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() - std::forward<Self>(self).getRHS().values();
    }

    template<Vector V, Scalar U>
    [[nodiscard, gnu::always_inline]] auto operator-(V&& v, U&& x) noexcept requires(!DeviceObj<V>) {
        return std::forward<V>(v) + (-std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard, gnu::always_inline]] auto operator-(U&& x, V&& v) noexcept requires(!DeviceObj<V>) {
        return std::forward<U>(x) + (-std::forward<V>(v));
    }

    template<Vector V1, Vector V2>
    [[nodiscard, gnu::always_inline]] auto operator-(V1&& v1, V2&& v2) noexcept requires(!DeviceObj<V1> && !DeviceObj<V2>) {
        return VectorExpr<ExprID::Sub, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}

#ifdef PHYSICA_MKL
    #include "MKL/Sub.h"
#endif
