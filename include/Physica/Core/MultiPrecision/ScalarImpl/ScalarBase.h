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
    template<class Derived> class ScalarBase;
    template<class T> class ScalarRef;
    template<class T> class ScalarPtr;
    template<class T, DiffMode Mode, int Order = 1> class Diff;

    template<class T>
    struct is_scalar : public std::is_base_of<ScalarBase<T>, T> {};

    template<class T>
    struct is_scalar<ScalarBase<T>> {
        constexpr static bool value = true;
    };

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

        using PlainScalar = typename Traits<Derived>::PlainScalar;
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

        using ValueType = typename std::conditional<isForwardDiff || !isDifferentiable, PlainScalar, PlainScalar*>::type;
        using GradType1 = typename std::conditional<Order == 1, ValueType, Diff<PlainScalar, Mode, Order - 1>>::type;
    public:
        using GradType = typename std::conditional<isDifferentiable, GradType1, PlainStruct<void>>::type;
        template<int GradOrder>
        using GradRtnTy = typename std::conditional<Order == GradOrder, PlainScalar, Diff<PlainScalar, Mode, Order - GradOrder>>::type;
    public:
        /* Operators */
        __host__ __device__ inline bool operator>(float s) const noexcept;
        __host__ __device__ inline bool operator<(float s) const noexcept;
        __host__ __device__ inline bool operator>(double s) const noexcept;
        __host__ __device__ inline bool operator<(double s) const noexcept;
        template<class T>
        __host__ __device__ inline bool operator>(const ScalarBase<T>& s) const noexcept;
        template<class T>
        __host__ __device__ inline bool operator<(const ScalarBase<T>& s) const noexcept;
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
        [[nodiscard]] __host__ __device__  const PlainScalar& getValue() const noexcept {
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
    template<class T>
    __host__ __device__ inline bool ScalarBase<Derived>::operator>(const ScalarBase<T>& s) const noexcept {
        static_assert(!isComplex && !T::isComplex, "[Error]: Comparison between complex scalars is invalid");
        return getValue() > s.getValue();
    }

    template<class Derived>
    template<class T>
    __host__ __device__ inline bool ScalarBase<Derived>::operator<(const ScalarBase<T>& s) const noexcept {
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

    template<class ScalarType>
    ScalarType relativeError(const ScalarType& scalar1, const ScalarType& scalar2) {
        static_assert(Core::is_scalar<ScalarType>::value && !ScalarType::isComplex && !ScalarType::isDifferentiable, "[Error]: Invalid template param");
        const auto& s1 = scalar1.getDerived();
        const auto& s2 = scalar2.getDerived();
        const ScalarType min = std::numeric_limits<ScalarType>::min();
        const bool useAbsCompare = (abs(s1) < min) || (abs(s2) < min);
        const ScalarType delta = s1 - s2;
        const ScalarType error = useAbsCompare ? abs(delta) : abs(delta / (s1 + s2) * ScalarType(2));
        return error;
    }

    template<class ScalarType>
    bool scalarNear(const ScalarBase<ScalarType>& s1,
                    const ScalarBase<ScalarType>& s2,
                    double precision) {
        assert(precision > 0);
        constexpr bool isDifferentiable = ScalarType::isDifferentiable;
        if constexpr (ScalarType::isComplex) {
            using PlainRealType = typename ScalarType::PlainScalar::RealType;
            const ScalarType diff = s1.getDerived() - s2.getDerived();
            const bool isValueNear = scalarNear(abs(diff.getValue()), PlainRealType(0), precision);
            if constexpr (isDifferentiable)
                return isValueNear && scalarNear(abs(diff.getGrad()), PlainRealType(0), precision);
            else
                return isValueNear;
        }
        else {
            using PlainScalar = typename ScalarType::PlainScalar;
            const bool isValueNear = relativeError(s1.getValue(), s2.getValue()) < PlainScalar(precision);
            if constexpr (ScalarType::isDifferentiable)
                return isValueNear && scalarNear(s1.getGrad(), s2.getGrad(), precision);
            else
                return isValueNear;
        }
    }

    template<class ScalarType>
    inline ScalarType operator+(const ScalarBase<ScalarType>& s) {
        return ScalarType(s.getDerived());
    }

    template<class ScalarType1, class ScalarType2>
    __host__ __device__ inline void operator+=(ScalarBase<ScalarType1>& s1, const ScalarBase<ScalarType2>& s2) {
        s1.getDerived() = s1.getDerived() + s2.getDerived();
    }

    template<class ScalarType1, class ScalarType2>
    __host__ __device__ inline void operator-=(ScalarBase<ScalarType1>& s1, const ScalarBase<ScalarType2>& s2) {
        s1.getDerived() = s1.getDerived() - s2.getDerived();
    }

    template<class ScalarType1, class ScalarType2>
    __host__ __device__ inline void operator*=(ScalarBase<ScalarType1>& s1, const ScalarBase<ScalarType2>& s2) {
        s1.getDerived() = s1.getDerived() * s2.getDerived();
    }

    template<class ScalarType1, class ScalarType2>
    __host__ __device__ inline void operator/=(ScalarBase<ScalarType1>& s1, const ScalarBase<ScalarType2>& s2) {
        s1.getDerived() = s1.getDerived() / s2.getDerived();
    }

    template<class ScalarType>
    inline void operator^=(ScalarBase<ScalarType>& s1, const ScalarBase<ScalarType>& s2) { s1.getDerived() = s1.getDerived() ^ s2.getDerived(); }

    template<class ScalarType>
    inline void operator<<=(ScalarBase<ScalarType>& s, int bits) { s.getDerived() = s.getDerived() << bits; }

    template<class ScalarType>
    inline void operator>>=(ScalarBase<ScalarType>& s, int bits) { s.getDerived() = s.getDerived() >> bits; }

    template<class ScalarType>
    __host__ __device__ inline bool operator>=(const ScalarBase<ScalarType>& s1, const ScalarBase<ScalarType>& s2) {
        return !(s1.getDerived() < s2.getDerived());
    }

    template<class ScalarType>
    __host__ __device__ inline bool operator<=(const ScalarBase<ScalarType>& s1, const ScalarBase<ScalarType>& s2) {
        return !(s1.getDerived() > s2.getDerived());
    }

    template<class ScalarType>
    __host__ __device__ inline bool operator!=(const ScalarBase<ScalarType>& s1, const ScalarBase<ScalarType>& s2) {
        return !(s1.getDerived() == s2.getDerived());
    }

    template<class ScalarType>
    [[nodiscard]] bool absCompare(const ScalarBase<ScalarType>& s1, const ScalarBase<ScalarType>& s2) {
        if constexpr (ScalarType::isComplex)
            return s1.squaredNorm() >= s2.squaredNorm();
        else
            return abs(s1.getDerived()) >= abs(s2.getDerived());
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
    template<class ScalarType>
    inline void swap(Physica::Core::ScalarBase<ScalarType>& __restrict s1, Physica::Core::ScalarBase<ScalarType>& __restrict s2) noexcept {
        s1.getDerived().swap(s2.getDerived());
    }
}
