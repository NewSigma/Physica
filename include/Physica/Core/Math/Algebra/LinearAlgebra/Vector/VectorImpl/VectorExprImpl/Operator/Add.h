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
    template<Vector V, Scalar U>
    class VectorExpr<ExprID::Add, V, U>
            : public BinaryVectorExpr<ExprID::Add, V, U> {
        using Base = BinaryVectorExpr<ExprID::Add, V, U>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs, size_t count) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const;

        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
        [[nodiscard]] __host__ __device__ consteval static bool lowerToFMA() noexcept;
    };

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Add, V, U>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept -> CoDiff<T> {
        if constexpr (lowerToFMA()) {
            T a = *(lhs.getLHS());
            T b;
            if constexpr (Scalar<decltype(lhs.getRHS())>)
                b = lhs.getRHS();
            else
                b = *(lhs.getRHS());
            T c = rhs;
            return fma(a, b, c);
        }
        else
            return *lhs + rhs;
    }

    template<Vector V, Scalar U>
    template<int Size>
    auto VectorExpr<ExprID::Add, V, U>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept -> SIMD<T, Size> {
        if constexpr (lowerToFMA()) {
            SIMD<T, Size> a = lhs.getLHS().template load<Size>();
            SIMD<T, Size> b;
            if constexpr (Scalar<decltype(lhs.getRHS())>)
                b = SIMD<T, Size>(lhs.getRHS());
            else
                b = lhs.getRHS().template load<Size>();
            auto c = SIMD<T, Size>(rhs);
            return fma(a, b, c);
        }
        else
            return lhs.template load<Size>() + SIMD<T, Size>(rhs);
    }

    template<Vector V, Scalar U>
    template<int Size>
    auto VectorExpr<ExprID::Add, V, U>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs, size_t count) noexcept -> SIMD<T, Size> {
        if constexpr (lowerToFMA()) {
            SIMD<T, Size> a = lhs.getLHS().template load<Size>(count);
            SIMD<T, Size> b;
            if constexpr (Scalar<decltype(lhs.getRHS())>)
                b = SIMD<T, Size>(lhs.getRHS());
            else
                b = lhs.getRHS().template load<Size>(count);
            auto c = SIMD<T, Size>(rhs, count);
            return fma(a, b, c);
        }
        else
            return lhs.template load<Size>(count) + SIMD<T, Size>(rhs, count);
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Add, V, U>::calc(size_t index) const -> CoDiff<T> {
        if constexpr (lowerToFMA()) {
            T a = getLHS().getLHS().calc(index);
            T b;
            if constexpr (Scalar<decltype(getLHS().getRHS())>)
                b = getLHS().getRHS();
            else
                b = getLHS().getRHS().calc(index);
            T c = getRHS();
            return fma(a, b, c);
        }
        else
            return getLHS().calc(index) + getRHS();
    }

    template<Vector V, Scalar U>
    void VectorExpr<ExprID::Add, V, U>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& g = grad.values();
        if constexpr (ReverseDiff<V>)
            Base::getLHS().reverse(g);
        if constexpr (ReverseDiff<U>)
            Base::getRHS().reverse(g.sum());
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Add, V, U>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() + std::forward<Self>(self).getRHS().value();
    }

    template<Vector V, Scalar U>
    consteval bool VectorExpr<ExprID::Add, V, U>::lowerToFMA() noexcept {
        if constexpr (instanceof_xt<VectorExpr, V>) {
            using V1 = std::remove_cvref_t<V>;
            using T1 = V1::ScalarType;
            using T2 = std::remove_cvref_t<U>;
            return (V1::getExprID() == ExprID::Mul) && std::same_as<T1, T2>;
        }
        return false;
    }

    template<Vector V1, Vector V2>
    class VectorExpr<ExprID::Add, V1, V2>
            : public BinaryVectorExpr<ExprID::Add, V1, V2> {
        using Base = BinaryVectorExpr<ExprID::Add, V1, V2>;
    public:
        using Base::isComplex;
        using Base::isReverseDiff;
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

        [[nodiscard]] CoDiff<T> calc(size_t index) const;

        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
        [[nodiscard]] __host__ __device__ consteval static bool lowerToFMA() noexcept;
    };

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Add, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept -> CoDiff<T> {
        if constexpr (lowerToFMA()) {
            T a = *(lhs.getLHS());
            T b;
            if constexpr (Scalar<decltype(lhs.getRHS())>)
                b = lhs.getRHS();
            else
                b = *(lhs.getRHS());
            T c = *rhs;
            return fma(a, b, c);
        }
        else
            return *lhs + *rhs;
    }

    template<Vector V1, Vector V2>
    template<int Size>
    auto VectorExpr<ExprID::Add, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept -> SIMD<T, Size> {
        if constexpr (lowerToFMA()) {
            SIMD<T, Size> a = lhs.getLHS().template load<Size>();
            SIMD<T, Size> b;
            if constexpr (Scalar<decltype(lhs.getRHS())>)
                b = SIMD<T, Size>(lhs.getRHS());
            else
                b = lhs.getRHS().template load<Size>();
            SIMD<T, Size> c = rhs.template load<Size>();
            return fma(a, b, c);
        }
        else
            return lhs.template load<Size>() + rhs.template load<Size>();
    }

    template<Vector V1, Vector V2>
    template<int Size>
    auto VectorExpr<ExprID::Add, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs, size_t count) noexcept -> SIMD<T, Size> {
        if constexpr (lowerToFMA()) {
            SIMD<T, Size> a = lhs.getLHS().template load<Size>(count);
            SIMD<T, Size> b;
            if constexpr (Scalar<decltype(lhs.getRHS())>)
                b = SIMD<T, Size>(lhs.getRHS());
            else
                b = lhs.getRHS().template load<Size>(count);
            SIMD<T, Size> c = rhs.template load<Size>(count);
            return fma(a, b, c);
        }
        else
            return lhs.template load<Size>(count) + rhs.template load<Size>(count);
    }

    template<Vector V1, Vector V2>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Add, V1, V2>::assign(Vector auto&& v) const {
        constexpr bool FastAssign1 = Traits<std::remove_cvref_t<V1>>::FastAssign;
        constexpr bool FastAssign2 = Traits<std::remove_cvref_t<V2>>::FastAssign;
        if constexpr (FastAssign1) {
            getLHS().template assign<P>(v);
            v += getRHS();
        }
        else if constexpr (FastAssign2) {
            getRHS().template assign<P>(v);
            v += getLHS();
        }
        else {
            using V = std::remove_cvref<decltype(v)>::type;
            constexpr size_t Size = std::max(Base::getSizeAtCompile(), v.getSizeAtCompile());
            constexpr size_t Critical = 256;
            constexpr bool UseMKL1 = Internal::EnableMKL<V1, V>::value;
            constexpr bool UseMKL2 = Internal::EnableMKL<V2, V>::value;
            constexpr bool UseMKL3 = Size == Dynamic || Size > Critical;
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
    auto VectorExpr<ExprID::Add, V1, V2>::calc(size_t index) const -> CoDiff<T> {
        if constexpr (lowerToFMA()) {
            T a = getLHS().getLHS().calc(index);
            T b;
            if constexpr (Scalar<decltype(getLHS().getRHS())>)
                b = getLHS().getRHS();
            else
                b = getLHS().getRHS().calc(index);
            T c = getRHS().calc(index);
            return fma(a, b, c);
        }
        else
            return getLHS().calc(index) + getRHS().calc(index);
    }

    template<Vector V1, Vector V2>
    void VectorExpr<ExprID::Add, V1, V2>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff());
        using U = decltype(grad);
        if constexpr (Scalar<U>) {
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(grad);
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(grad);
        }
        else {
            static_assert(Vector<U>);
            const auto& g = grad.values();
            assert(g.getLength() == Base::getLength());
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(g);
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(g);
        }
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Add, V1, V2>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() + std::forward<Self>(self).getRHS().values();
    }

    template<Vector V1, Vector V2>
    consteval bool VectorExpr<ExprID::Add, V1, V2>::lowerToFMA() noexcept {
        if constexpr (instanceof_xt<VectorExpr, V1>) {
            using T1 = std::remove_cvref_t<V1>::ScalarType;
            using T2 = std::remove_cvref_t<V2>::ScalarType;
            return (std::remove_cvref_t<V1>::getExprID() == ExprID::Mul) && std::same_as<T1, T2>;
        }
        return false;
    }

    template<Vector V, Scalar U>
    [[nodiscard, gnu::always_inline]] auto operator+(V&& v, U&& x) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Add, V&&, U&&>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard, gnu::always_inline]] auto operator+(U&& x, V&& v) noexcept requires(!DeviceObj<V>) {
        return std::forward<V>(v) + std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard, gnu::always_inline]] auto operator+(V1&& v1, V2&& v2) noexcept requires(!DeviceObj<V1> && !DeviceObj<V2>) {
        if constexpr (!canonicalized(v1, v2))
            return std::forward<V2>(v2) + std::forward<V1>(v1);
        else
            return VectorExpr<ExprID::Add, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}

#ifdef PHYSICA_MKL
    #include "MKL/Add.h"
#endif
