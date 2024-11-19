/*
 * Copyright 2023-2024 Weibo He.
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
#include <Physica/CRTPBase.h>
#include <Physica/PlainStruct.h>
#include <Physica/Core/Exception/NoImplException.h>

namespace Physica::Core {
    enum class DiffMode {
        Forward,
        Reverse
    };
    template<class ScalarType> class ScalarRef;
    template<class ScalarType> class ScalarPtr;
    template<class ScalarType, DiffMode Mode, int Order = 1> class Diff;

    template<class Derived>
    class ScalarBase : public CRTPBase<ScalarBase<Derived>> {
        using Base = CRTPBase<ScalarBase<Derived>>;
    public:
        constexpr static ScalarOption Option = Traits<Derived>::Option;
        constexpr static int Order = Traits<Derived>::Order;
        constexpr static bool isComplex = Traits<Derived>::isComplex;
        constexpr static bool isDifferentiable = Traits<Derived>::isDifferentiable;
        constexpr static bool isForwardDiff = Traits<Derived>::isForwardDiff;
        constexpr static bool isReverseDiff = Traits<Derived>::isReverseDiff;
        constexpr static DiffMode Mode = isForwardDiff ? DiffMode::Forward : DiffMode::Reverse;

        using ValueType = typename Traits<Derived>::ValueType;
        using ScalarType = typename Traits<Derived>::ScalarType;
        using PtrTy = typename Traits<Derived>::PtrTy;
        using ConstPtrTy = typename Traits<Derived>::ConstPtrTy;
        using RefTy = typename Traits<Derived>::RefTy;
        using ConstRefTy = typename Traits<Derived>::ConstRefTy;
        using RealType = typename Traits<Derived>::RealType;
        using ComplexType = typename Traits<Derived>::ComplexType;
        using MachineType = typename Traits<Derived>::MachineType;
        using device_obj_type = Derived;
    private:
        constexpr static bool isConsistent1 = isDifferentiable && (isForwardDiff != isReverseDiff) && (Order > 0);
        constexpr static bool isConsistent2 = !isDifferentiable && !isForwardDiff && !isReverseDiff && (Order == 0);
        static_assert(isConsistent1 || isConsistent2, "[Error]: DiffMode is not self consistent");
        static_assert(std::is_same<Derived, ScalarType>::value, "[Error]: Inconsistence type between traits and inherit class");

        using ValuePtr = typename std::conditional<isForwardDiff || !isDifferentiable, ValueType, ValueType*>::type;
        using GradType1 = typename std::conditional<Order == 1, ValuePtr, Diff<ValueType, Mode, Order - 1>>::type;
    public:
        using GradType = typename std::conditional<isDifferentiable, GradType1, PlainStruct<void>>::type;
        template<int GradOrder>
        using GradRtnTy = typename std::conditional<Order == GradOrder, ValueType, Diff<ValueType, Mode, Order - GradOrder>>::type;
    public:
        /* Operators */
        __host__ __device__ inline bool operator>(float s) const noexcept;
        __host__ __device__ inline bool operator<(float s) const noexcept;
        __host__ __device__ inline bool operator>(double s) const noexcept;
        __host__ __device__ inline bool operator<(double s) const noexcept;
        template<Scalar T>
        __host__ __device__ inline bool operator>(const T& s) const noexcept;
        template<Scalar T>
        __host__ __device__ inline bool operator<(const T& s) const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ RealType real() const;
        [[nodiscard]] __host__ __device__ RealType imag() const;
        [[nodiscard]] __host__ __device__ ScalarType conjugate() const {
            if constexpr (isComplex)
                return this->getDerived().conjugate();
            else
                return real();
        }

        [[nodiscard]] __host__ __device__ ScalarType unit() const {
            if constexpr (isComplex)
                return this->getDerived().unit();
            else
                return ScalarType(real().isNegative() ? -1 : 1);
        }

        [[nodiscard]] __host__ __device__ RealType norm() const { return sqrt(squaredNorm()); }

        [[nodiscard]] __host__ __device__ RealType squaredNorm() const {
            if constexpr (isComplex)
                return this->getDerived().squaredNorm();
            else
                return square(this->getDerived());
        }
        /* SIMD support */
        __host__ __device__ constexpr static size_t size() { return 1; }

        __host__ __device__ Derived& load(const ScalarType* p) {
            this->getDerived() = *p;
            return this->getDerived();
        }

        __host__ __device__ void store(ScalarType* p) const { *p = this->getDerived().getValue(); }

        __host__ __device__ Derived& load_partial(const ScalarType* p, int n) {
            if (n)
                load(p);
            return this->getDerived();
        }

        __host__ __device__ void store_partial(ScalarType* p, int n) const {
            if (n)
                store(p);
        }

        void insert([[maybe_unused]] int index, ScalarType value) { this->getDerived() = ScalarType(value); }
        [[nodiscard]] ScalarType sum() const { return this->getDerived(); }
        /* Auto differential support */
        [[nodiscard]] __host__ __device__  const ValueType& getValue() const noexcept {
            if constexpr (isDifferentiable)
                return this->getDerived().getValue();
            else
                return this->getDerived();
        }

        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__  const GradRtnTy<GradOrder>& getGrad() const noexcept {
            if constexpr (isDifferentiable)
                return this->getDerived().getGrad();
            else
                noImpl();
        }
        /* Static Members */
        static inline bool matchSign(const ScalarType& s1, const ScalarType& s2) {
            return (s1.isPositive() && s2.isPositive()) || (s1.isNegative() && s2.isNegative());
        }
    };

    template<class Derived>
    __host__ __device__ inline bool ScalarBase<Derived>::operator>(float s) const noexcept {
        static_assert(!isComplex, "[Error]: Comparison between complex scalars is invalid");
        return float(getValue()) > s;
    }

    template<class Derived>
    __host__ __device__ inline bool ScalarBase<Derived>::operator<(float s) const noexcept {
        static_assert(!isComplex, "[Error]: Comparison between complex scalars is invalid");
        return float(getValue()) < s;
    }

    template<class Derived>
    __host__ __device__ inline bool ScalarBase<Derived>::operator>(double s) const noexcept {
        static_assert(!isComplex, "[Error]: Comparison between complex scalars is invalid");
        return double(getValue()) > s;
    }

    template<class Derived>
    __host__ __device__ inline bool ScalarBase<Derived>::operator<(double s) const noexcept {
        static_assert(!isComplex, "[Error]: Comparison between complex scalars is invalid");
        return double(getValue()) < s;
    }
    /**
     * TODO: Move to \class Diff once we merge forward and reverse pass
     */
    template<class Derived>
    template<Scalar T>
    __host__ __device__ inline bool ScalarBase<Derived>::operator>(const T& s) const noexcept {
        static_assert(!isComplex && !T::isComplex, "[Error]: Comparison between complex scalars is invalid");
        return getValue() > s.getValue();
    }

    template<class Derived>
    template<Scalar T>
    __host__ __device__ inline bool ScalarBase<Derived>::operator<(const T& s) const noexcept {
        static_assert(!isComplex && !T::isComplex, "[Error]: Comparison between complex scalars is invalid");
        return getValue() < s.getValue();
    }

    template<class Derived>
    __host__ __device__ typename ScalarBase<Derived>::RealType ScalarBase<Derived>::real() const {
        if constexpr (isDifferentiable)
            return RealType(getValue().real(), getGrad().real());
        else if constexpr (isComplex)
            return this->getDerived().real();
        else
            return this->getDerived();
    }

    template<class Derived>
    __host__ __device__ typename ScalarBase<Derived>::RealType ScalarBase<Derived>::imag() const {
        if constexpr (isDifferentiable)
            return RealType(getValue().imag(), getGrad().imag());
        else if constexpr (isComplex)
            return this->getDerived().imag();
        else
            return RealType(0);
    }

    template<Scalar T>
    T relativeError(const T& x, const T& y) {
        static_assert(!T::isComplex && !T::isDifferentiable, "[Error]: Invalid template param");
        const T min = std::numeric_limits<T>::min();
        const bool useAbsCompare = (abs(x) < min) || (abs(y) < min);
        const T delta = x - y;
        const T error = useAbsCompare ? abs(delta) : abs(delta / (x + y) * T(2));
        return error;
    }

    template<Scalar T>
    bool scalarNear(const T& x, const T& y, double precision) {
        assert(precision > 0);
        constexpr bool isDifferentiable = T::isDifferentiable;
        if constexpr (T::isComplex) {
            using PlainRealType = typename T::ValueType::RealType;
            const T diff = x - y;
            const bool isValueNear = scalarNear(abs(diff.getValue()), PlainRealType(0), precision);
            if constexpr (isDifferentiable)
                return isValueNear && scalarNear(abs(diff.getGrad()), PlainRealType(0), precision);
            else
                return isValueNear;
        }
        else {
            using ValueType = typename T::ValueType;
            const bool isValueNear = relativeError(x.getValue(), y.getValue()) < ValueType(precision);
            if constexpr (T::isDifferentiable)
                return isValueNear && scalarNear(x.getGrad(), y.getGrad(), precision);
            else
                return isValueNear;
        }
    }

    template<Scalar T>
    inline T operator+(const T& x) {
        return T(x);
    }

    template<Scalar T1, Scalar T2>
    __host__ __device__ inline void operator+=(T1& x, const T2& y) {
        x = x + y;
    }

    template<Scalar T1, Scalar T2>
    __host__ __device__ inline void operator-=(T1& x, const T2& y) {
        x = x - y;
    }

    template<Scalar T1, Scalar T2>
    __host__ __device__ inline void operator*=(T1& x, const T2& y) {
        x = x * y;
    }

    template<Scalar T1, Scalar T2>
    __host__ __device__ inline void operator/=(T1& x, const T2& y) {
        x = x / y;
    }

    template<Scalar T>
    inline void operator^=(T& x, const T& y) { x = x ^ y; }

    template<Scalar T>
    inline void operator<<=(T& x, int bits) { x = x << bits; }

    template<Scalar T>
    inline void operator>>=(T& x, int bits) { x = x >> bits; }

    template<Scalar T>
    __host__ __device__ inline bool operator>=(const T& x, const T& y) {
        return !(x < y);
    }

    template<Scalar T>
    __host__ __device__ inline bool operator<=(const T& x, const T& y) {
        return !(x > y);
    }

    template<Scalar T>
    __host__ __device__ inline bool operator!=(const T& x, const T& y) {
        return !(x == y);
    }

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
    class Traits<Core::ScalarBase<T>> {
    public:
        using Derived = T;
    };
}

namespace std {
    template<Physica::Core::Scalar T>
    inline void swap(T& __restrict x, T& __restrict y) noexcept {
        x.swap(y);
    }
}
