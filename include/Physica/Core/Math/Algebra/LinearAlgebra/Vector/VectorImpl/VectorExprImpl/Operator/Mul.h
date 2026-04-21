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
    class VectorExpr<ExprID::Mul, V, U>
            : public BinaryVectorExpr<ExprID::Mul, V, U> {
        using Base = BinaryVectorExpr<ExprID::Mul, V, U>;
    public:
        using Base::isComplex;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        using Base::operator*;
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept;
        [[nodiscard]] static CoDiff<T> operator()(const Scalar auto& lhs, std::random_access_iterator auto rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs, size_t count) noexcept;
        [[nodiscard]] auto operator*(Scalar auto x) const noexcept;
        [[nodiscard]] auto operator-() const& noexcept;
        [[nodiscard]] auto operator-() && noexcept;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const;
        void assign_mkl(Vector auto& v) const noexcept;

        template<ExecutePolicy P = Sequential>
        void assign_add(Vector auto&& v) const;
        void assign_add_mkl(Vector auto&& v) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_add_base(Vector auto&& v) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] T sum() const { return getLHS().sum() * getRHS(); }

        [[nodiscard]] auto values(this auto&& self) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    private:
        template<ExecutePolicy P>
        void assign_fma_for(Vector auto& __restrict v) const __restrict noexcept;
        template<Vector Target, size_t Length>
        void assign_fma_simd(Target& v) const noexcept;
    };

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Mul, V, U>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept -> CoDiff<T> {
        return (*lhs) * rhs;
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Mul, V, U>::operator()(const Scalar auto& lhs, std::random_access_iterator auto rhs) noexcept -> CoDiff<T> {
        return lhs * (*rhs);
    }

    template<Vector V, Scalar U>
    template<int Size>
    auto VectorExpr<ExprID::Mul, V, U>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept -> SIMD<T, Size> {
        using V1 = std::remove_cvref_t<V>;
        using U1 = std::remove_cvref_t<U>;
        if constexpr (!V1::isComplex() && U1::isComplex()) {
            auto pack = lhs.template load<Size>();
            return SIMD<T, Size>::asComplex(SIMD<Tr, Size * 2>(pack * rhs.real(), pack * rhs.imag()).gatherRealImag());
        }
        else if constexpr (V1::isComplex() && !U1::isComplex()) {
            auto pack = lhs.template load<Size>();
            return SIMD<T, Size>::asComplex(pack.asReal() * SIMD<Tr, Size * 2>(rhs));
        }
        else
            return lhs.template load<Size>() * SIMD<T, Size>(rhs);
    }

    template<Vector V, Scalar U>
    template<int Size>
    auto VectorExpr<ExprID::Mul, V, U>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs, size_t count) noexcept -> SIMD<T, Size> {
        using V1 = std::remove_cvref_t<V>;
        using U1 = std::remove_cvref_t<U>;
        if constexpr (!V1::isComplex() && U1::isComplex()) {
            auto pack = lhs.template load<Size>(count);
            return SIMD<T, Size>::asComplex(SIMD<Tr, Size * 2>(pack * rhs.real(), pack * rhs.imag()).gatherRealImag());
        }
        else if constexpr (V1::isComplex() && !U1::isComplex()) {
            auto pack = lhs.template load<Size>(count);
            return SIMD<T, Size>::asComplex(pack.asReal() * SIMD<Tr, Size * 2>(rhs));
        }
        else
            return lhs.template load<Size>(count) * SIMD<T, Size>(rhs);
    }

    // FIXME: we cannot use explicit object parameter, clang 22 rejects because ambiguous overloads
    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Mul, V, U>::operator*(Scalar auto x) const noexcept {
        return getLHS() * (getRHS() * x);
    }

    // FIXME: we cannot use explicit object parameter because of regression; seems clang 22 has problem on overload?
    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Mul, V, U>::operator-() const& noexcept {
        return getLHS() * (-getRHS());
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Mul, V, U>::operator-() && noexcept {
        return std::move(getLHS()) * (-getRHS());
    }

    template<Vector V, Scalar U>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Mul, V, U>::assign(Vector auto&& v) const {
        if constexpr (v.isFastAssign()) {
            getLHS().template assign<P>(v);
            v *= getRHS();
        }
        else
            Base::assign(v);
    }

    template<Vector V, Scalar U>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Mul, V, U>::assign_add(Vector auto&& v) const {
        using V1 = std::remove_cvref_t<decltype(v)>;
        constexpr size_t Size = std::max(Base::getSizeAtCompile(), v.getSizeAtCompile());
        constexpr bool SmallVector = 0 < Size && Size <= 128;
        if constexpr (Internal::EnableMKL<V, V1>::value && !SmallVector) {
            if (Base::getLength() <= 128)
                assign_add_base(v);
            else
                assign_add_mkl(v);
        }
        else
            assign_add_base(v);
    }

    template<Vector V, Scalar U>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Mul, V, U>::assign_add_base(Vector auto&& v) const noexcept {
        using Source = std::remove_cvref<V>::type;
        using Target = std::remove_cvref<decltype(v)>::type;
        using T1 = Source::ScalarType;
        using T2 = Target::ScalarType;
        constexpr bool LowerToFMA = std::same_as<T1, T2>;
        if constexpr (LowerToFMA) {
            v.assert_assign(Base::getDerived());
            if constexpr (Internal::EnableSIMD<Source, Target>::value && !isReverseDiff()) {
                constexpr size_t Length = std::max(Base::getSizeAtCompile(), v.getSizeAtCompile());
                assign_fma_simd<Target, Length>(v);
            }
            else
                assign_fma_for<P>(v);
        }
        else
            Base::template assign_add<P>(v);
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Mul, V, U>::calc(size_t index) const -> CoDiff<T> {
        return getLHS().calc(index) * getRHS();
    }

    template<Vector V, Scalar U>
    void VectorExpr<ExprID::Mul, V, U>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& g = grad.values();
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        if constexpr (ReverseDiff<V>)
            lhs.reverse(rhs.value() * g);
        if constexpr (ReverseDiff<U>)
            rhs.reverse(lhs.values() * g);
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Mul, V, U>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().value();
    }

    template<Vector V, Scalar U>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Mul, V, U>::assign_fma_for(Vector auto& __restrict v) const __restrict noexcept {
        parallel_for<P>([&, this](size_t i) {
            if constexpr (isReverseDiff())
                v[i] = fma(getLHS().calc_value(i), Tv(getRHS().value()), v[i]);
            else
                v[i] = fma(getLHS().calc(i), T(getRHS().value()), v[i]);
        }, Base::getLength(), 0).wait();
    }

    template<Vector V, Scalar U>
    template<Vector Target, size_t Length>
    void VectorExpr<ExprID::Mul, V, U>::assign_fma_simd(Target& v) const noexcept {
        using Pack = BestPacket<typename Target::ScalarType, Length>::Type;
        constexpr size_t Size = Pack::size();
        const auto rhs = Pack(getRHS());
        const T b = getRHS();
        auto it = zip(v.view(), getLHS().view()).begin();
        if constexpr (Length != Dynamic) {
            constexpr size_t to = Length / Size * Size;
            for (size_t i = 0; i < to; i += Size) {
                auto [it_v, it_lhs] = it + i;
                it_v.store(fma(it_lhs.template load<Size>(), rhs, it_v.template load<Size>()));
            }

            for (size_t i = Length - Length % Size; i < Length; ++i) {
                auto [it_v, it_lhs] = it + i;
                *it_v = fma(*it_lhs, b, *it_v);
            }
        }
        else {
            const size_t length = v.getLength();
            size_t i = 0;
            for (; i < length / Size * Size; i += Size) {
                auto [it_v, it_lhs] = it + i;
                it_v.store(fma(it_lhs.template load<Size>(), rhs, it_v.template load<Size>()));
            }

            for (; i < length; ++i) {
                auto [it_v, it_lhs] = it + i;
                *it_v = fma(*it_lhs, b, *it_v);
            }
        }
    }
    ////////////////////////////////////////////////////////////////////
    template<Vector V1, Vector V2>
    class VectorExpr<ExprID::Mul, V1, V2>
            : public BinaryVectorExpr<ExprID::Mul, V1, V2> {
        using This = VectorExpr<ExprID::Mul, V1, V2>;
        using Base = BinaryVectorExpr<ExprID::Mul, V1, V2>;
    public:
        using Base::isComplex;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tr;
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

        template<ExecutePolicy P = Sequential>
        void assign_add(Vector auto&& v) const;
        template<ExecutePolicy P = Sequential>
        void assign_add_base(Vector auto&& v) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        using Base::reverse;
        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] T sum(this auto&&) noexcept;

        [[nodiscard]] auto values(this auto&& self) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    private:
        template<ExecutePolicy P>
        void assign_fma_for(Vector auto& __restrict v) const __restrict noexcept;
        template<Vector Target, size_t Size>
        void assign_fma_simd(Target& v) const noexcept;
    };

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Mul, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept -> CoDiff<T> {
        return (*lhs) * (*rhs);
    }

    template<Vector V1, Vector V2>
    template<int Size>
    auto VectorExpr<ExprID::Mul, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept -> SIMD<T, Size> {
        using LHS = std::remove_cvref_t<V1>;
        using RHS = std::remove_cvref_t<V2>;
        auto pack1 = lhs.template load<Size>();
        auto pack2 = rhs.template load<Size>();
        if constexpr (LHS::isComplex() && !RHS::isComplex())
            return SIMD<T, Size>::asComplex(pack1.asReal() * SIMD<Tr, 2 * Size>(pack2, pack2).scatterRealImag());
        else
            return pack1 * pack2;
    }

    template<Vector V1, Vector V2>
    template<int Size>
    auto VectorExpr<ExprID::Mul, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs, size_t count) noexcept -> SIMD<T, Size> {
        using LHS = std::remove_cvref_t<V1>;
        using RHS = std::remove_cvref_t<V2>;
        auto pack1 = lhs.template load<Size>(count);
        auto pack2 = rhs.template load<Size>(count);
        if constexpr (LHS::isComplex() && !RHS::isComplex())
            return SIMD<T, Size>::asComplex(pack1.asReal() * SIMD<Tr, 2 * Size>(pack2, pack2).scatterRealImag());
        else
            return pack1 * pack2;
    }

    template<Vector V1, Vector V2>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Mul, V1, V2>::assign(Vector auto&& v) const {
        Base::template assign_base<P>(v);
    }

    template<Vector V1, Vector V2>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Mul, V1, V2>::assign_add(Vector auto&& v) const {
        assign_add_base<P>(std::forward<decltype(v)>(v));
    }

    template<Vector V1, Vector V2>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Mul, V1, V2>::assign_add_base(Vector auto&& v) const noexcept {
        using Target = std::remove_cvref_t<decltype(v)>;
        using T1 = std::remove_cvref_t<V1>::ScalarType;
        using T2 = std::remove_cvref_t<V2>::ScalarType;
        using U = Target::ScalarType;
        constexpr bool LowerToFMA = std::same_as<T1, T2> && std::same_as<T, U>;
        if constexpr (LowerToFMA) {
            v.assert_assign(Base::getDerived());
            if constexpr (Internal::EnableSIMD<This, Target>::value && !isReverseDiff()) {
                constexpr size_t Size = std::max(Base::getSizeAtCompile(), v.getSizeAtCompile());
                assign_fma_simd<Target, Size>(v);
            }
            else
                assign_fma_for<P>(v);
        }
        else
            Base::template assign_add<P>(v);
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Mul, V1, V2>::calc(size_t index) const -> CoDiff<T> {
        return getLHS().calc(index) * getRHS().calc(index);
    }

    template<Vector V1, Vector V2>
    void VectorExpr<ExprID::Mul, V1, V2>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff());
        if constexpr (Scalar<decltype(grad)>) {
            if constexpr (ReverseDiff<V1>)
                getLHS().reverse(getRHS().values() * grad);
            if constexpr (ReverseDiff<V2>)
                getRHS().reverse(getLHS().values() * grad);
        }
        else {
            static_assert(Vector<decltype(grad)>);
            const auto& g = grad.values();
            if constexpr (ReverseDiff<V1>)
                getLHS().reverse(hadamard(getRHS().values(), g));
            if constexpr (ReverseDiff<V2>)
                getRHS().reverse(hadamard(getLHS().values(), g));
        }
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Mul, V1, V2>::sum(this auto&& self) noexcept -> T {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS() * std::forward<Self>(self).getRHS();
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Mul, V1, V2>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return hadamard(std::forward<Self>(self).getLHS().values(), std::forward<Self>(self).getRHS().values());
    }

    template<Vector V1, Vector V2>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Mul, V1, V2>::assign_fma_for(Vector auto& __restrict v) const __restrict noexcept {
        parallel_for<P>([&, this](size_t i) {
            if constexpr (isReverseDiff())
                v[i] = fma(getLHS().calc_value(i), Tv(getRHS().calc_value(i)), v[i]);
            else
                v[i] = fma(getLHS().calc(i), T(getRHS().calc(i)), v[i]);
        }, Base::getLength(), 0).wait();
    }

    template<Vector V1, Vector V2>
    template<Vector Target, size_t Length>
    void VectorExpr<ExprID::Mul, V1, V2>::assign_fma_simd(Target& v) const noexcept {
        using Pack = BestPacket<typename Target::ScalarType, Length>::Type;
        constexpr size_t Size = Pack::size();
        auto it = zip(v.view(), getLHS().view(), getRHS().view()).begin();
        if constexpr (Length != Dynamic) {
            constexpr size_t to = Length / Size * Size;
            for (size_t i = 0; i < to; i += Size) {
                auto [it_v, it_lhs, it_rhs] = it + i;
                it_v.store(fma(it_lhs.template load<Size>(), it_rhs.template load<Size>(), it_v.template load<Size>()));
            }

            for (size_t i = Length - Length % Size; i < Length; ++i) {
                auto [it_v, it_lhs, it_rhs] = it + i;
                *it_v = fma(*it_lhs, *it_rhs, *it_v);
            }
        }
        else {
            const size_t length = v.getLength();
            size_t i = 0;
            for (; i < length / Size * Size; i += Size) {
                auto [it_v, it_lhs, it_rhs] = it + i;
                it_v.store(fma(it_lhs.template load<Size>(), it_rhs.template load<Size>(), it_v.template load<Size>()));
            }

            for (; i < length; ++i) {
                auto [it_v, it_lhs, it_rhs] = it + i;
                *it_v = fma(*it_lhs, *it_rhs, *it_v);
            }
        }
    }
    ////////////////////////////////////////////////////////////////////
    template<Vector V, Scalar U>
    [[nodiscard, gnu::always_inline]] auto operator*(V&& v, U&& x) noexcept requires(!DeviceObj<V>) {
        using RtnTy = VectorExpr<ExprID::Mul, V&&, U&&>;
        if constexpr (instanceof_xt<VectorExpr, V>) {
            using RHS = Traits<V>::RHS;
            if constexpr (v.getExprID() == ExprID::Mul && Scalar<RHS>)
                return v.getLHS() * (v.getRHS() * x);
            else
                return RtnTy(std::forward<V>(v), std::forward<U>(x));
        }
        else
            return RtnTy(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard, gnu::always_inline]] auto operator*(U&& x, V&& v) noexcept requires(!DeviceObj<V>) {
        return std::forward<V>(v) * std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard, gnu::always_inline]] auto hadamard(V1&& v1, V2&& v2) noexcept requires(!DeviceObj<V1> && !DeviceObj<V2>) {
        using RtnTy = VectorExpr<ExprID::Mul, V1&&, V2&&>;
        if constexpr (!canonicalized(v1, v2))
            return hadamard(std::forward<V2>(v2), std::forward<V1>(v1));
        else if constexpr (instanceof_xt<VectorExpr, V1>) {
            using RHS1 = Traits<V1>::RHS;
            if constexpr (v1.getExprID() == ExprID::Mul && Scalar<RHS1>) // if we can lower V1
                return hadamard(v2, v1.getLHS()) * v1.getRHS();
            else if constexpr (instanceof_xt<VectorExpr, V2>) { // if not, see if we can lower V2
                using RHS2 = Traits<V2>::RHS;
                if constexpr (v2.getExprID() == ExprID::Mul && Scalar<RHS2>)
                    return hadamard(v1, v2.getLHS()) * v2.getRHS();
                else
                    return RtnTy(std::forward<V1>(v1), std::forward<V2>(v2));
            }
            else // there is nothing we can do here
                return RtnTy(std::forward<V1>(v1), std::forward<V2>(v2));
        }
        else {
            static_assert(!instanceof_xt<VectorExpr, V2>, "[Error]: Canonicalization failed");
            return RtnTy(std::forward<V1>(v1), std::forward<V2>(v2));
        }
    }
}

#include "MKL/Mul.h"
