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

#include <iomanip>
#include "Physica/Utils/Template/CRTPBase.h"

namespace Physica::Core {
    enum class DiffMode {
        Forward,
        Reverse
    };
    template<class Derived> class ScalarBase;
    template<class T, DiffMode Mode, unsigned int Order> class Differentiable;

    template<class T>
    struct is_scalar : public std::is_base_of<ScalarBase<T>, T> {};

    template<class T>
    struct is_scalar<ScalarBase<T>> {
        constexpr static bool value = true;
    };

    template<class Derived>
    class ScalarBase : public Utils::CRTPBase<Derived> {
    public:
        using ScalarType = typename Traits<Derived>::ScalarType;
        using TrivialType = typename Traits<Derived>::TrivialType;
        using RealType = typename Traits<Derived>::RealType;
        using ComplexType = typename Traits<Derived>::ComplexType;
        using PlainScalar = typename Traits<Derived>::PlainScalar;
        constexpr static ScalarOption Option = Traits<Derived>::Option;
        constexpr static bool isComplex = Traits<Derived>::isComplex;
        constexpr static bool isDifferentiable = Traits<Derived>::isDifferentiable;
        constexpr static bool isForwardDiff = Traits<Derived>::isForwardDiff;
        constexpr static bool isReverseDiff = Traits<Derived>::isReverseDiff;
        constexpr static unsigned int Order = Traits<Derived>::Order;
    private:
        constexpr static bool isConsistent1 = isDifferentiable && (isForwardDiff != isReverseDiff) && (Order > 0);
        constexpr static bool isConsistent2 = !isDifferentiable && !isForwardDiff && !isReverseDiff && (Order == 0);
        static_assert(isConsistent1 || isConsistent2, "[Error]: DiffMode is not self consistent");
        static_assert(std::is_same<Derived, ScalarType>::value, "[Error]: Inconsistence type between traits and inherit class");
    public:
        /* Operators */
        template<class T>
        __host__ __device__ inline bool operator>(const ScalarBase<T>& s) const {
            static_assert(!isComplex && !T::isComplex, "[Error]: Comparison between complex scalars is invalid");
            return getValue() > s.getValue();
        }

        template<class T>
        __host__ __device__ inline bool operator<(const ScalarBase<T>& s) const {
            static_assert(!isComplex && !T::isComplex, "[Error]: Comparison between complex scalars is invalid");
            return getValue() < s.getValue();
        }
        /* Getters */
        [[nodiscard]] __host__ __device__ const RealType& getReal() const {
            if constexpr (isComplex)
                return this->getDerived().getReal();
            else
                return this->getDerived();
        }

        [[nodiscard]] __host__ __device__ RealType getImag() const {
            if constexpr (isComplex)
                return this->getDerived().getImag();
            else
                return Derived(0);
        }

        [[nodiscard]] __host__ __device__ ScalarType conjugate() const {
            if constexpr (isComplex)
                return this->getDerived().conjugate();
            else
                return getReal();
        }

        [[nodiscard]] __host__ __device__ ScalarType unit() const {
            if constexpr (isComplex)
                return this->getDerived().unit();
            else
                return ScalarType(getReal().isNegative() ? -1 : 1);
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

        Derived& load(const ScalarType* p) {
            this->getDerived() = *p;
            return this->getDerived();
        }

        void store(ScalarType* p) const { *p = this->getDerived().getValue(); }

        Derived& load_partial(int n, const ScalarType* p) {
            if (n)
                load(p);
            return this->getDerived();
        }

        void store_partial(int n, ScalarType* p) const {
            if (n)
                store(p);
        }

        void insert([[maybe_unused]] int index, ScalarType value) { this->getDerived() = ScalarType(value); }
        [[nodiscard]] ScalarType horizontal_add() const { return this->getDerived(); }
        /* Auto differential support */
        [[nodiscard]] const PlainScalar& getValue() const noexcept {
            if constexpr (isDifferentiable)
                return this->getDerived().getValue();
            else
                return this->getDerived();
        }

        [[nodiscard]] const PlainScalar& getGrad() const noexcept {
            return this->getDerived().getGrad();
        }
    };

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
            const bool isValueNear = relativeError(s1.getValue().getReal(), s2.getValue().getReal()) < PlainScalar(precision);
            if constexpr (ScalarType::isDifferentiable)
                return isValueNear && relativeError(s1.getGrad().getReal(), s2.getGrad().getReal()) < PlainScalar(precision);
            else
                return isValueNear;
        }
    }

    template<ScalarOption Option>
    std::ostream& operator<<(std::ostream& os, const Scalar<Option>& s) {
        constexpr int FloatPrec = 7;
        constexpr int DoublePrec = 16;
        const auto lastPrec = os.precision();
        return os << std::setprecision(Option == Float ? FloatPrec : DoublePrec)
                  << double(s)
                  << std::setprecision(lastPrec);
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

    template<class ScalarType>
    inline void swap(ScalarBase<ScalarType>& __restrict s1, ScalarBase<ScalarType>& __restrict s2) noexcept {
        s1.getDerived().swap(s2.getDerived());
    }
}
