/*
 * Copyright 2023-2026 Weibo He.
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

#include <cassert>
#include <cmath>
#include <type_traits>
#include "Physica/CRTPBase.h"
#include "Physica/Core/Scalar/Scalar.h"

namespace Physica {
    namespace Internal {
        template<class T, int GradOrder>
        class GradTypeHelper {
            static_assert(Traits<T>::Order >= GradOrder, "[Error]: Required order is too high");
            static_assert(GradOrder > 1, "[Error]: Bad GradOrder");
        public:
            using Type = GradTypeHelper<typename T::GradType, GradOrder - 1>::Type;
        };

        template<class T>
        class GradTypeHelper<T, 1> {
            constexpr static int Order = Traits<T>::Order;
            constexpr static DiffMode Mode = Traits<T>::isForwardDiff ? DiffMode::Forward : DiffMode::Reverse;
            static_assert(Order >= 1, "[Error]: Required order is too high");
            using ValueType = Traits<T>::ValueType;
        public:
            using Type = std::conditional<Order == 1, ValueType, Diff<ValueType, Mode, Order - 1>>::type;
        };

        template<class T>
        class GradTypeHelper<T, 0> {
        public:
            using Type = void;
        };
    }

    template<class Derived>
    class ScalarBase : public CRTPBase<ScalarBase<Derived>> {
        using This = ScalarBase<Derived>;
        using Base = CRTPBase<This>;
        using TraitsType = Traits<Derived>;
    public:
        constexpr static FloatPrec Prec = TraitsType::Prec;
        constexpr static int Order = TraitsType::Order;
        constexpr static DiffMode Mode = TraitsType::isForwardDiff ? DiffMode::Forward : DiffMode::Reverse;

        using ScalarType = TraitsType::ScalarType;
        using ValueType = TraitsType::ValueType;
        using GradType = Internal::GradTypeHelper<ScalarType, Order == 0 ? 0 : 1>::Type;
        using PtrTy = TraitsType::PtrTy;
        using ConstPtrTy = TraitsType::ConstPtrTy;
        using RefTy = TraitsType::RefTy;
        using ConstRefTy = TraitsType::ConstRefTy;
        using RealType = TraitsType::RealType;
        using ComplexType = TraitsType::ComplexType;
        using MachineType = TraitsType::MachineType;
        using device_obj_type = Derived;

        template<int GradOrder>
        using GradWithOrder = Internal::GradTypeHelper<ScalarType, Order == 0 ? 0 : GradOrder>::Type;
        template<int MaskOrder>
        using GradMaskType = std::conditional<MaskOrder == 0, ValueType, Diff<ValueType, Mode, MaskOrder>>::type;
    public:
        constexpr ~ScalarBase() = default;
        /* Operators */
        __host__ __device__ Derived& operator=(const Scalar auto& x);
        __host__ __device__ void operator+=(this auto&, const Scalar auto& x) noexcept;
        __host__ __device__ void operator-=(this auto&, const Scalar auto& x) noexcept;
        __host__ __device__ void operator*=(this auto&, const Scalar auto& x) noexcept;
        __host__ __device__ void operator/=(this auto&, const Scalar auto& x) noexcept;
        [[nodiscard]] __host__ __device__ auto operator<=>(this const auto& x, float y) noexcept;
        [[nodiscard]] __host__ __device__ auto operator<=>(this const auto& x, double y) noexcept;
        [[nodiscard]] __host__ __device__ auto operator<=>(this const auto& x, const Scalar auto& y) noexcept;
        /* Operations */
        [[nodiscard]] __host__ __device__ ScalarType conjugate() const noexcept;
        [[nodiscard]] __host__ __device__ RealType norm() const noexcept;
        [[nodiscard]] __host__ __device__ auto squaredNorm() const noexcept;
        [[nodiscard]] __host__ __device__ ScalarType sum() const noexcept;

        __host__ __device__ void toNextMean(size_t lastNumSample, ScalarType sample) noexcept;
        __host__ __device__ void toNextVariance(ScalarType& __restrict mean, size_t lastNumSample, ScalarType sample) noexcept;
        template<class R>
        void random_uniform() noexcept;
        template<class R>
        void random_normal() noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ decltype(auto) real(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ decltype(auto) imag(this auto&&) noexcept;

        [[nodiscard]] __host__ __device__ auto* real_ptr(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto* value_ptr(this auto&&) noexcept;
        template<int GradOrder>
        [[nodiscard]] __host__ __device__ auto* grad_ptr(this auto&&) noexcept;

        [[nodiscard]] __host__ __device__ auto& value(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ auto& grad(this auto&&) noexcept;
        template<int MaskOrder>
        [[nodiscard]] __host__ __device__ decltype(auto) grad_mask(this auto&&) noexcept;

        [[nodiscard]] __host__ __device__ auto isZero() const noexcept;
        [[nodiscard]] __host__ __device__ auto isSubNormal() const noexcept;
        [[nodiscard]] __host__ __device__ auto isPositive() const noexcept;
        [[nodiscard]] __host__ __device__ auto isNegative() const noexcept;
        [[nodiscard]] __host__ __device__ auto isFinite() const noexcept;
        /* Static Members */
        [[nodiscard]] __host__ __device__ consteval static bool isComplex() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isDiffable() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isForwardDiff() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isReverseDiff() noexcept;
        static bool matchSign(const ScalarType& s1, const ScalarType& s2);
        template<Scalar Src>
        consteval static void static_assert_assign() noexcept;
        consteval static void checkComplexCompare() noexcept;
    protected:
        constexpr ScalarBase() = default;
        constexpr ScalarBase(const This&) = default;
        constexpr ScalarBase(This&&) = default;
        /* Operators */
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
    private:
        static_assert(isDiffable() == (Order > 0), "[Error]: DiffMode is not self consistent");
        static_assert(std::is_same<Derived, ScalarType>::value || instanceof<Derived, ScalarRef>, "[Error]: Inconsistent type between traits and inherit class");
    };

    template<class Derived>
    __host__ __device__ Derived& ScalarBase<Derived>::operator=(const Scalar auto& x) {
        static_assert_assign<decltype(x)>();
        return Base::getDerived() = ScalarType(x);
    }

    template<class Derived>
    __host__ __device__ void ScalarBase<Derived>::operator+=(this auto& self, const Scalar auto& x) noexcept {
        static_assert_assign<decltype(x)>();
        self = self + x;
    }

    template<class Derived>
    __host__ __device__ void ScalarBase<Derived>::operator-=(this auto& self, const Scalar auto& x) noexcept {
        static_assert_assign<decltype(x)>();
        self = self - x;
    }

    template<class Derived>
    __host__ __device__ void ScalarBase<Derived>::operator*=(this auto& self, const Scalar auto& x) noexcept {
        static_assert_assign<decltype(x)>();
        self = self * x;
    }

    template<class Derived>
    __host__ __device__ void ScalarBase<Derived>::operator/=(this auto& self, const Scalar auto& x) noexcept {
        static_assert_assign<decltype(x)>();
        self = self / x;
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::operator<=>(this const auto& x, float y) noexcept {
        checkComplexCompare();
        return float(x.value()) <=> y;
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::operator<=>(this const auto& x, double y) noexcept {
        checkComplexCompare();
        return double(x.value()) <=> y;
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::operator<=>(this const auto& x, const Scalar auto& y) noexcept {
        using X = std::remove_cvref<decltype(x)>::type;
        using Y = std::remove_cvref<decltype(y)>::type;
        static_assert(X::isDiffable() || Y::isDiffable() || !std::same_as<X, Y>, "[Error]: This function handles type casts only");
        using Z = Internal::BinaryScalarOpRtnTy<typename X::ValueType, typename Y::ValueType>::Type;
        return Z(x.value()) <=> Z(y.value());
    }

    template<class Derived>
    __host__ __device__ decltype(auto) ScalarBase<Derived>::real(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isDiffable())
            return RealType(self.value().real(), self.grad().real());
        else if constexpr (isComplex())
            return std::forward<Self>(self).real();
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    __host__ __device__ decltype(auto) ScalarBase<Derived>::imag(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isDiffable())
            return RealType(self.value().imag(), self.grad().imag());
        else if constexpr (isComplex())
            return std::forward<Self>(self).imag();
        else
            return RealType(0);
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::conjugate() const noexcept -> ScalarType {
        const auto& x = this->getDerived();
        if constexpr (isDiffable())
            return ScalarType(x.value().conjugate(), x.grad().conjugate());
        else if constexpr (isComplex())
            return ScalarType(x.real(), -x.imag());
        else
            return x.real();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::norm() const noexcept -> RealType {
        return sqrt(squaredNorm());
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::squaredNorm() const noexcept {
        if constexpr (isComplex() || isReverseDiff())
            return this->getDerived().squaredNorm();
        else
            return square(this->getDerived());
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::sum() const noexcept -> ScalarType {
        return this->getDerived();
    }
    /**
     * Algorithms of toNext*() are also known as welford accumulators
     */
    template<class Derived>
    __host__ __device__ void ScalarBase<Derived>::toNextMean(size_t lastNumSample, ScalarType sample) noexcept {
        auto& mean = Base::getDerived();
        if (lastNumSample == 0) [[unlikely]] {
            mean = sample;
            return;
        }
        mean = fma(ScalarType(lastNumSample), mean, sample) / ScalarType(lastNumSample + 1);
    }

    template<class Derived>
    __host__ __device__ void ScalarBase<Derived>::toNextVariance(ScalarType& __restrict mean, size_t lastNumSample, ScalarType sample) noexcept {
        assert(this != &mean && "[Error]: Var and mean cannot overlap");
        auto& var = Base::getDerived();
        if (lastNumSample == 0) [[unlikely]] {
            mean = sample;
            var = 0;
            return;
        }
        const auto factor = reciprocal(ScalarType(lastNumSample + 1));
        var = fma(square(mean - sample), factor, var) * (ScalarType(lastNumSample) * factor);
        mean.toNextMean(lastNumSample, sample);
    }

    template<class Derived>
    template<class R>
    void ScalarBase<Derived>::random_uniform() noexcept {
        Base::getDerived() = Derived::template random_uniform<R>();
    }

    template<class Derived>
    template<class R>
    void ScalarBase<Derived>::random_normal() noexcept {
        Base::getDerived() = Derived::template random_normal<R>();
    }

    template<class Derived>
    __host__ __device__ auto* ScalarBase<Derived>::real_ptr(this auto&& self) noexcept {
        if constexpr (isComplex())
            return self.real_ptr();
        else
            return &self;
    }

    template<class Derived>
    __host__ __device__ auto* ScalarBase<Derived>::value_ptr(this auto&& self) noexcept {
        if constexpr (isDiffable())
            return self.value_ptr();
        else
            return &self;
    }

    template<class Derived>
    template<int GradOrder>
    __host__ __device__ auto* ScalarBase<Derived>::grad_ptr(this auto&& self) noexcept {
        static_assert(isDiffable(), "[Error]: Cannot take grad ptr of a undiffable scalar");
        return self.template grad_ptr<GradOrder>();
    }

    template<class Derived>
    __host__ __device__ auto& ScalarBase<Derived>::value(this auto&& self) noexcept {
        return *self.value_ptr();
    }

    template<class Derived>
    template<int GradOrder>
    __host__ __device__  auto& ScalarBase<Derived>::grad(this auto&& self) noexcept {
        return *self.template grad_ptr<GradOrder>();
    }

    template<class Derived>
    template<int MaskOrder>
    __host__ __device__ decltype(auto) ScalarBase<Derived>::grad_mask(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isDiffable())
            return std::forward<Self>(self).template grad_mask<MaskOrder>();
        else
            return std::forward<Self>(self); // Provides a syntactic sugar for non-diffable scalar
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::isZero() const noexcept {
        if constexpr (isDiffable())
            return Base::getDerived().value().isZero();
        else
            return Base::getDerived().isZero();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::isSubNormal() const noexcept {
        if constexpr (isDiffable())
            return Base::getDerived().value().isSubNormal();
        else
            return Base::getDerived().isSubNormal();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::isPositive() const noexcept {
        checkComplexCompare();
        if constexpr (isDiffable())
            return Base::getDerived().value().isPositive();
        else
            return Base::getDerived().isPositive();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::isNegative() const noexcept {
        checkComplexCompare();
        if constexpr (isDiffable())
            return Base::getDerived().value().isNegative();
        else
            return Base::getDerived().isNegative();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::isFinite() const noexcept {
        return Base::getDerived().isFinite();
    }

    template<class Derived>
    __host__ __device__ consteval bool ScalarBase<Derived>::isComplex() noexcept {
        return TraitsType::isComplex;
    }

    template<class Derived>
    __host__ __device__ consteval bool ScalarBase<Derived>::isForwardDiff() noexcept {
        return TraitsType::isForwardDiff;
    }

    template<class Derived>
    __host__ __device__ consteval bool ScalarBase<Derived>::isReverseDiff() noexcept {
        return TraitsType::isReverseDiff;
    }

    template<class Derived>
    __host__ __device__ consteval bool ScalarBase<Derived>::isDiffable() noexcept {
        return isForwardDiff() || isReverseDiff();
    }

    template<class Derived>
    bool ScalarBase<Derived>::matchSign(const ScalarType& s1, const ScalarType& s2) {
        return (s1.isPositive() && s2.isPositive()) || (s1.isNegative() && s2.isNegative());
    }

    template<class Derived>
    template<Scalar Src>
    consteval void ScalarBase<Derived>::static_assert_assign() noexcept {
        using U = std::remove_cvref_t<Src>;
        static_assert(!ReverseDiff<Src>, "[Error]: Assign reverse diffable scalar to another scalar breaks compute graph");
        static_assert(isDiffable() || !U::isDiffable(), "[Error]: Assign diffable scalar to normal scalar discards grad");
        static_assert(isComplex() || !U::isComplex(), "[Error]: Assign complex to real discards imag");
    }

    template<class Derived>
    consteval void ScalarBase<Derived>::checkComplexCompare() noexcept {
        static_assert(!isComplex(), "[Error]: Comparison between complex scalars is invalid");
    }

    template<Scalar T>
    T relativeError(T x, T y) noexcept {
        static_assert(!T::isComplex() && !T::isDiffable(), "[Error]: Invalid template param");
        T error = abs(x - y);
        if (x.isSubNormal() || y.isSubNormal())
            return error;

        T mean = (abs(x) + abs(y)) * 0.5;
        return error / mean;
    }

    template<Scalar T>
    [[clang::no_sanitize("numerical")]] bool scalarNear(const T& x, const T& y, double precision) noexcept {
        assert(precision > 0);
        if constexpr (T::isDiffable())
            return scalarNear(x.value(), y.value(), precision) && scalarNear(x.grad(), y.grad(), precision);
        else {
            using Tv = T::ValueType;
            if constexpr (T::isComplex()) {
                using Trv = Tv::RealType;
                Trv diff = std::hypot(x.value().real().toMachine() - y.value().real().toMachine(), x.value().imag().toMachine() - y.value().imag().toMachine());
                return scalarNear(diff, Trv(0), precision);
            }
            else
                return relativeError(x.value(), y.value()) < Tv(precision);
        }
    }

    template<Scalar T>
    bool scalarNear(const T& x, const T& y, uint64_t ulp) noexcept {
        if constexpr (T::isDiffable())
            return scalarNear(x.value(), y.value(), ulp) && scalarNear(x.grad(), y.grad(), ulp);
        else {
            if constexpr (T::isComplex())
                return scalarNear(x.real(), y.real(), ulp)
                    && scalarNear(x.imag(), y.imag(), ulp);
            else
                return calcULPDiff(x, y) <= ulp;
        }
    }

    template<Scalar T>
    [[nodiscard]] T&& operator+(T&& x) noexcept {
        return std::forward<T>(x);
    }

    template<Scalar T>
    void operator^=(T& x, const T& y) { x = x ^ y; }

    template<Scalar T>
    void operator<<=(T& x, int bits) { x = x << bits; }

    template<Scalar T>
    void operator>>=(T& x, int bits) { x = x >> bits; }

    template<Scalar T>
    [[nodiscard]] bool absCompare(const T& x, const T& y) {
        if constexpr (T::isComplex())
            return x.squaredNorm() >= y.squaredNorm();
        else
            return abs(x) >= abs(y);
    }

    template<Scalar T>
    void swap(T& x, T& y) noexcept {
        x.swap(y);
    }
}

namespace Physica {
    template<class T>
    class Traits<ScalarBase<T>> {
    public:
        using Derived = T;
    };
}
