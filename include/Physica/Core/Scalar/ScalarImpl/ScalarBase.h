/*
 * Copyright 2023-2025 Weibo He.
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
        constexpr static bool isComplex = TraitsType::isComplex;
        constexpr static bool isForwardDiff = TraitsType::isForwardDiff;
        constexpr static bool isReverseDiff = TraitsType::isReverseDiff;
        constexpr static bool isDiffable = isForwardDiff || isReverseDiff;
        constexpr static DiffMode Mode = isForwardDiff ? DiffMode::Forward : DiffMode::Reverse;

        using ValueType = TraitsType::ValueType;
        using ScalarType = TraitsType::ScalarType;
        using PtrTy = TraitsType::PtrTy;
        using ConstPtrTy = TraitsType::ConstPtrTy;
        using RefTy = TraitsType::RefTy;
        using ConstRefTy = TraitsType::ConstRefTy;
        using RealType = TraitsType::RealType;
        using ComplexType = TraitsType::ComplexType;
        using MachineType = TraitsType::MachineType;
        using device_obj_type = Derived;

        using GradType = Internal::GradTypeHelper<ScalarType, isDiffable ? 1 : 0>::Type;
        using MKL_Complex = std::conditional<Prec == Float32, MKL_Complex8,
                                                              typename std::conditional<Prec == Float64, MKL_Complex16, void>::type>::type;
    private:
        constexpr static bool isScalarRef = Internal::IsScalarRef<Derived>::value;
        static_assert(isDiffable == (Order > 0), "[Error]: DiffMode is not self consistent");
        static_assert(std::is_same<Derived, ScalarType>::value || isScalarRef, "[Error]: Inconsistent type between traits and inherit class");

        template<int GradOrder>
        using GradRtnTy = std::conditional<!isDiffable || Order == GradOrder, ValueType, Diff<ValueType, Mode, Order - GradOrder>>::type;
    public:
        constexpr ~ScalarBase() = default;
        /* Operators */
        template<Scalar U>
        __host__ __device__ Derived& operator=(const U& x) requires(isReverseDiff || !ReverseDiff<U>);
        template<Scalar U>
        __host__ __device__ void operator+=(const U& x) requires(isReverseDiff || !ReverseDiff<U>);
        template<Scalar U>
        __host__ __device__ void operator-=(const U& x) requires(isReverseDiff || !ReverseDiff<U>);
        template<Scalar U>
        __host__ __device__ void operator*=(const U& x) requires(isReverseDiff || !ReverseDiff<U>);
        template<Scalar U>
        __host__ __device__ void operator/=(const U& x) requires(isReverseDiff || !ReverseDiff<U>);
        __host__ __device__ bool operator>(float s) const noexcept;
        __host__ __device__ bool operator<(float s) const noexcept;
        __host__ __device__ bool operator>(double s) const noexcept;
        __host__ __device__ bool operator<(double s) const noexcept;
        template<Scalar T>
        __host__ __device__ bool operator>(const T& x) const noexcept;
        template<Scalar T>
        __host__ __device__ bool operator<(const T& x) const noexcept;
        __host__ __device__ bool operator>=(const Scalar auto& x) const noexcept;
        __host__ __device__ bool operator<=(const Scalar auto& x) const noexcept;
        __host__ __device__ bool operator!=(const Scalar auto& x) const noexcept;
        /* Operations */
        [[nodiscard]] ScalarType calc() const;

        template<int MaskOrder>
        [[nodiscard]] __host__ __device__ auto mask() const noexcept;

        __host__ __device__ Derived& load(ConstPtrTy p);
        __host__ __device__ void store(PtrTy p) const;
        __host__ __device__ Derived& load_partial(ConstPtrTy p, int n);
        __host__ __device__ void store_partial(PtrTy p, int n) const;
        void insert(int index, ScalarType value);

        [[nodiscard]] __host__ __device__ RealType real() const noexcept;
        [[nodiscard]] __host__ __device__ RealType imag() const noexcept;
        [[nodiscard]] __host__ __device__ ScalarType conjugate() const noexcept;
        [[nodiscard]] __host__ __device__ RealType norm() const noexcept;
        [[nodiscard]] __host__ __device__ auto squaredNorm() const noexcept;
        [[nodiscard]] __host__ __device__ ScalarType sum() const noexcept;
        /* Getters */
        __host__ __device__ constexpr static size_t size() { return 1; }
        [[nodiscard]] __host__ __device__ auto value_ptr() noexcept;
        [[nodiscard]] __host__ __device__ const auto value_ptr() const noexcept;
        template<int GradOrder>
        [[nodiscard]] __host__ __device__ auto grad_ptr() noexcept;
        template<int GradOrder>
        [[nodiscard]] __host__ __device__ const auto grad_ptr() const noexcept;

        [[nodiscard]] __host__ __device__ ValueType& value() noexcept;
        [[nodiscard]] __host__ __device__ const ValueType& value() const noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ GradRtnTy<GradOrder>& grad() noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ const GradRtnTy<GradOrder>& grad() const noexcept;
        [[nodiscard]] __host__ __device__ auto isZero() const noexcept;
        [[nodiscard]] __host__ __device__ auto isPositive() const noexcept;
        [[nodiscard]] __host__ __device__ auto isNegative() const noexcept;
        /* Static Members */
        static bool matchSign(const ScalarType& s1, const ScalarType& s2);
        template<Scalar U>
        consteval static void checkAssign();
        consteval static void checkComplexCompare();
    protected:
        constexpr ScalarBase() = default;
        constexpr ScalarBase(const This&) = default;
        constexpr ScalarBase(This&&) = default;
        /* Operators */
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
    };

    template<class Derived>
    template<Scalar U>
    __host__ __device__ Derived& ScalarBase<Derived>::operator=(const U& x) requires(isReverseDiff || !ReverseDiff<U>) {
        checkAssign<U>();
        return Base::getDerived() = ScalarType(x);
    }

    template<class Derived>
    template<Scalar U>
    __host__ __device__ void ScalarBase<Derived>::operator+=(const U& x) requires(isReverseDiff || !ReverseDiff<U>) {
        checkAssign<U>();
        auto& y = Base::getDerived();
        if constexpr (ReverseDiff<U>)
            y += x.value();
        else if constexpr (isReverseDiff)
            y.value() += x;
        else
            y = y + x;
    }

    template<class Derived>
    template<Scalar U>
    __host__ __device__ void ScalarBase<Derived>::operator-=(const U& x) requires(isReverseDiff || !ReverseDiff<U>) {
        checkAssign<U>();
        auto& y = Base::getDerived();
        if constexpr (ReverseDiff<U>)
            y -= x.value();
        else if constexpr (isReverseDiff)
            y.value() -= x;
        else
            y = y - x;
    }

    template<class Derived>
    template<Scalar U>
    __host__ __device__ void ScalarBase<Derived>::operator*=(const U& x) requires(isReverseDiff || !ReverseDiff<U>) {
        checkAssign<U>();
        auto& y = Base::getDerived();
        if constexpr (ReverseDiff<U>)
            y *= x.value();
        else if constexpr (isReverseDiff)
            y.value() *= x;
        else
            y = y * x;
    }

    template<class Derived>
    template<Scalar U>
    __host__ __device__ void ScalarBase<Derived>::operator/=(const U& x) requires(isReverseDiff || !ReverseDiff<U>) {
        checkAssign<U>();
        auto& y = Base::getDerived();
        if constexpr (ReverseDiff<U>)
            y /= x.value();
        else if constexpr (isReverseDiff)
            y.value() /= x;
        else
            y = y / x;
    }

    template<class Derived>
    __host__ __device__ bool ScalarBase<Derived>::operator>(float s) const noexcept {
        checkComplexCompare();
        return float(value()) > s;
    }

    template<class Derived>
    __host__ __device__ bool ScalarBase<Derived>::operator<(float s) const noexcept {
        checkComplexCompare();
        return float(value()) < s;
    }

    template<class Derived>
    __host__ __device__ bool ScalarBase<Derived>::operator>(double s) const noexcept {
        checkComplexCompare();
        return double(value()) > s;
    }

    template<class Derived>
    __host__ __device__ bool ScalarBase<Derived>::operator<(double s) const noexcept {
        checkComplexCompare();
        return double(value()) < s;
    }
    /**
     * TODO: Move to \class Diff once we merge forward and reverse pass
     */
    template<class Derived>
    template<Scalar T>
    __host__ __device__ bool ScalarBase<Derived>::operator>(const T& x) const noexcept {
        static_assert(!isComplex && !T::isComplex, "[Error]: Comparison between complex scalars is invalid");
        return value() > x.value();
    }

    template<class Derived>
    template<Scalar T>
    __host__ __device__ bool ScalarBase<Derived>::operator<(const T& x) const noexcept {
        static_assert(!isComplex && !T::isComplex, "[Error]: Comparison between complex scalars is invalid");
        return value() < x.value();
    }

    template<class Derived>
    __host__ __device__ bool ScalarBase<Derived>::operator>=(const Scalar auto& x) const noexcept {
        return !operator<(x);
    }

    template<class Derived>
    __host__ __device__ bool ScalarBase<Derived>::operator<=(const Scalar auto& x) const noexcept {
        return !operator>(x);
    }

    template<class Derived>
    __host__ __device__ bool ScalarBase<Derived>::operator!=(const Scalar auto& x) const noexcept {
        return !operator==(x);
    }

    template<class Derived>
    auto ScalarBase<Derived>::calc() const -> ScalarType {
        return ScalarType(Base::getDerived());
    }

    template<class Derived>
    template<int MaskOrder>
    __host__ __device__ auto ScalarBase<Derived>::mask() const noexcept {
        if constexpr (isDiffable)
            return calc().template mask<MaskOrder>();
        else
            return Base::getDerived();
    }

    template<class Derived>
    __host__ __device__ Derived& ScalarBase<Derived>::load(ConstPtrTy p) {
        this->getDerived() = (*p).value();
        return this->getDerived();
    }

    template<class Derived>
    __host__ __device__ void ScalarBase<Derived>::store(PtrTy p) const {
        *p = this->getDerived().value();
    }

    template<class Derived>
    __host__ __device__ Derived& ScalarBase<Derived>::load_partial(ConstPtrTy p, int n) {
        if (n)
            load(p);
        return this->getDerived();
    }

    template<class Derived>
    __host__ __device__ void ScalarBase<Derived>::store_partial(PtrTy p, int n) const {
        if (n)
            store(p);
    }

    template<class Derived>
    void ScalarBase<Derived>::insert([[maybe_unused]] int index, ScalarType value) {
        assert(index == 0);
        this->getDerived() = ScalarType(value);
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::real() const noexcept -> RealType {
        const auto& x = this->getDerived();
        if constexpr (isDiffable)
            return RealType(x.value().real(), x.grad().real());
        else if constexpr (isComplex)
            return x.real();
        else
            return x;
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::imag() const noexcept -> RealType {
        const auto& x = this->getDerived();
        if constexpr (isDiffable)
            return RealType(x.value().imag(), x.grad().imag());
        else if constexpr (isComplex)
            return x.imag();
        else
            return RealType(0);
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::conjugate() const noexcept -> ScalarType {
        const auto& x = this->getDerived();
        if constexpr (isDiffable)
            return ScalarType(x.value().conjugate(), x.grad().conjugate());
        else if constexpr (isComplex)
            return ScalarType(real(), -imag());
        else
            return real();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::norm() const noexcept -> RealType {
        return sqrt(squaredNorm());
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::squaredNorm() const noexcept {
        if constexpr (isComplex || isReverseDiff)
            return this->getDerived().squaredNorm();
        else
            return square(this->getDerived());
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::sum() const noexcept -> ScalarType {
        return this->getDerived();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::value_ptr() noexcept {
        if constexpr (isDiffable)
            return Base::getDerived().value_ptr();
        else
            return &Base::getDerived();
    }

    template<class Derived>
    __host__ __device__ const auto ScalarBase<Derived>::value_ptr() const noexcept {
        return Base::getConstCastDerived().value_ptr();
    }

    template<class Derived>
    template<int GradOrder>
    __host__ __device__ auto ScalarBase<Derived>::grad_ptr() noexcept {
        static_assert(isDiffable, "[Error]: Cannot take grad ptr of a undiffable scalar");
        return Base::getDerived().template grad_ptr<GradOrder>();
    }

    template<class Derived>
    template<int GradOrder>
    __host__ __device__ const auto ScalarBase<Derived>::grad_ptr() const noexcept {
        return Base::getConstCastDerived().template grad_ptr<GradOrder>();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::value() noexcept -> ValueType& {
        return *value_ptr();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::value() const noexcept -> const ValueType& {
        return *value_ptr();
    }

    template<class Derived>
    template<int GradOrder>
    __host__ __device__  auto ScalarBase<Derived>::grad() noexcept -> GradRtnTy<GradOrder>& {
        return *grad_ptr<GradOrder>();
    }

    template<class Derived>
    template<int GradOrder>
    __host__ __device__  auto ScalarBase<Derived>::grad() const noexcept -> const GradRtnTy<GradOrder>& {
        return *grad_ptr<GradOrder>();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::isZero() const noexcept {
        if constexpr (isDiffable)
            return value().isZero();
        else
            return Base::getDerived().isZero();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::isPositive() const noexcept {
        checkComplexCompare();
        if constexpr (isDiffable)
            return value().isPositive();
        else
            return Base::getDerived().isPositive();
    }

    template<class Derived>
    __host__ __device__ auto ScalarBase<Derived>::isNegative() const noexcept {
        checkComplexCompare();
        if constexpr (isDiffable)
            return value().isNegative();
        else
            return Base::getDerived().isNegative();
    }

    template<class Derived>
    bool ScalarBase<Derived>::matchSign(const ScalarType& s1, const ScalarType& s2) {
        return (s1.isPositive() && s2.isPositive()) || (s1.isNegative() && s2.isNegative());
    }

    template<class Derived>
    template<Scalar U>
    consteval void ScalarBase<Derived>::checkAssign() {
        static_assert(isDiffable || !U::isDiffable, "[Error]: Assign diffable scalar to normal scalar discards grad");
        static_assert(isComplex || !U::isComplex, "[Error]: Assign complex to real discards imag");
    }

    template<class Derived>
    consteval void ScalarBase<Derived>::checkComplexCompare() {
        static_assert(!isComplex, "[Error]: Comparison between complex scalars is invalid");
    }

    template<Scalar T>
    T relativeError(const T& x, const T& y) noexcept {
        static_assert(!T::isComplex && !T::isDiffable, "[Error]: Invalid template param");
        const T min = std::numeric_limits<T>::min();
        T absX = abs(x);
        T absY = abs(y);
        T error = abs(x - y);
        bool useAbsCompare = (absX < min) || (absY < min);
        if (!useAbsCompare)
            error *= T(2) / (absX + absY);
        return error;
    }

    template<Scalar T>
    bool scalarNear(const T& x, const T& y, double precision) noexcept {
        assert(precision > 0);
        constexpr bool isDiffable = T::isDiffable;
        if constexpr (T::isComplex) {
            using PlainRealType = T::ValueType::RealType;
            const T diff = x - y;
            const bool isValueNear = scalarNear(abs(diff.value()), PlainRealType(0), precision);
            if constexpr (isDiffable)
                return isValueNear && scalarNear(abs(diff.grad()), PlainRealType(0), precision);
            else
                return isValueNear;
        }
        else {
            using ValueType = T::ValueType;
            const bool isValueNear = relativeError(x.value(), y.value()) < ValueType(precision);
            if constexpr (T::isDiffable)
                return isValueNear && scalarNear(x.grad(), y.grad(), precision);
            else
                return isValueNear;
        }
    }

    template<Scalar T>
    T&& operator+(T&& x) noexcept {
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
        if constexpr (T::isComplex)
            return x.squaredNorm() >= y.squaredNorm();
        else
            return abs(x) >= abs(y);
    }
}

namespace Physica {
    template<class T>
    class Traits<ScalarBase<T>> {
    public:
        using Derived = T;
    };
}

namespace std {
    template<Physica::Scalar T>
    void swap(T& __restrict x, T& __restrict y) noexcept {
        x.swap(y);
    }
}
