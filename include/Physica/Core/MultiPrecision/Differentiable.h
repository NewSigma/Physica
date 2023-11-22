/*
 * Copyright 2023 WeiBo He.
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

#include "Scalar.h"
#include "DifferentiableImpl/DiffTracer.h"

namespace Physica::Core {
    namespace Internal {
        template<class T, DiffMode M>
        class Traits<Differentiable<T, M>> {
            static_assert(!Traits<T>::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");
            using RealT = typename T::RealType;
            using ComplexT = typename T::ComplexType;
        public:
            using ScalarType = Differentiable<T, M>;
            using RealType = Differentiable<RealT, M>;
            using ComplexType = Differentiable<ComplexT, M>;
            using TrivialType = typename T::TrivialType;
            using PlainScalar = T;
            constexpr static ScalarOption Option = T::Option;
            constexpr static bool isComplex = T::isComplex;
            constexpr static bool isDifferentiable = true;
            constexpr static DiffMode Mode = M;
            /* SIMD */
            using BoolSIMDType = BoolSIMD<ScalarType, 1>;
        };

        template<class ScalarType, DiffMode Mode, class OtherScalar>
        class BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar> {
        public:
            using Type = Differentiable<typename BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type, Mode>;
        };

        template<class ScalarType, DiffMode Mode, class OtherScalar>
        class BinaryScalarOpReturnType<OtherScalar, Differentiable<ScalarType, Mode>> {
        public:
            using Type = typename BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type;
        };

        template<class ScalarType, DiffMode Mode, class OtherScalar, DiffMode OtherMode>
        class BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, Differentiable<OtherScalar, OtherMode>> {
            static_assert(Mode == OtherMode, "[Error]: DiffMode do not match");
        public:
            using Type = Differentiable<typename BinaryScalarOpReturnType<ScalarType, OtherScalar>::Type, Mode>;
        };
    }
    /**
     * \class Differentiable provides auto differential support for scalars
     */
    template<class ScalarType>
    class Differentiable<ScalarType, DiffMode::Forward> : public ScalarBase<Differentiable<ScalarType, DiffMode::Forward>> {
        using This = Differentiable<ScalarType, DiffMode::Forward>;
        using DiffTracerType = DiffTracer<ScalarType>;

        ScalarType value;
        ScalarType tangent;
    public:
        Differentiable() = default;
        Differentiable(double d) : This(ScalarType(d)) {}
        Differentiable(ScalarType value_);
        Differentiable(ScalarType value_, ScalarType tangent_);
        Differentiable(const Differentiable&) = default;
        Differentiable(Differentiable&&) noexcept = default;
        ~Differentiable() = default;
        /* Operators */
        Differentiable& operator=(Differentiable obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] explicit operator float() const { return float(value); }
        [[nodiscard]] explicit operator double() const { return double(value); }
        [[nodiscard]] inline bool operator==(const This& other) const;
        [[nodiscard]] inline Differentiable operator-() const;
        /* Operations */
        void swap(Differentiable& obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType& getValue() noexcept { return value; }
        [[nodiscard]] ScalarType& getTangent() noexcept { return tangent; }
        [[nodiscard]] const ScalarType& getValue() const noexcept { return value; }
        [[nodiscard]] const ScalarType& getTangent() const noexcept { return tangent; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return value.isZero(); }
        [[nodiscard]] __host__ __device__ bool isPositive() const { return value.isPositive(); }
        [[nodiscard]] __host__ __device__ bool isNegative() const { return value.isNegative(); }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] inline static Differentiable random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static Differentiable random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static Differentiable random_any(Distribution& dist, RandomGenerator& gen);
    };
    ////////////////////////////////////////////////////////////
    template<class ScalarType>
    class Differentiable<ScalarType, DiffMode::Reverse> : public ScalarBase<Differentiable<ScalarType, DiffMode::Reverse>> {
        using This = Differentiable<ScalarType, DiffMode::Reverse>;
        using DiffTracerType = DiffTracer<ScalarType>;

        ScalarType* __restrict pValue;
        ScalarType* __restrict pTangent;
    public:
        Differentiable() = default;
        Differentiable(double d) : This(ScalarType(d)) {}
        Differentiable(ScalarType value);
        Differentiable(ScalarType value, ScalarType tangent);
        Differentiable(ScalarType* pValue_, ScalarType* pTangent_);
        Differentiable(const Differentiable&) = default;
        Differentiable(Differentiable&&) noexcept = default;
        ~Differentiable() = default;
        /* Operators */
        Differentiable& operator=(const Differentiable&) = default;
        Differentiable& operator=(Differentiable&&) noexcept = default;
        [[nodiscard]] explicit operator float() const { return float(getValue()); }
        [[nodiscard]] explicit operator double() const { return double(getValue()); }
        [[nodiscard]] inline bool operator==(const This& other) const;
        [[nodiscard]] inline Differentiable operator-() const;
        /* Operations */
        void reverse();
        void reverse(Differentiable to);
        [[nodiscard]] Differentiable copy() const;
        inline void swap(Differentiable& obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType* value_ptr() const noexcept { return pValue; }
        [[nodiscard]] ScalarType* tangent_ptr() const noexcept { return pTangent; }
        [[nodiscard]] ScalarType& getValue() noexcept { return *pValue; }
        [[nodiscard]] ScalarType& getTangent() noexcept { return *pTangent; }
        [[nodiscard]] const ScalarType& getValue() const noexcept { return *pValue; }
        [[nodiscard]] const ScalarType& getTangent() const noexcept { return *pTangent; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return getValue().isZero(); }
        [[nodiscard]] __host__ __device__ bool isPositive() const { return getValue().isPositive(); }
        [[nodiscard]] __host__ __device__ bool isNegative() const { return getValue().isNegative(); }
        [[nodiscard]] inline ExpressionType getSource() const noexcept;
        /* Setters */
        void setValue(const ScalarType& x) { *pValue = x; }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] inline static Differentiable random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static Differentiable random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static Differentiable random_any(Distribution& dist, RandomGenerator& gen);
    };
    ////////////////////////////////////////////////////////////
    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator+(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator+(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2);

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator-(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator-(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2);

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator*(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator*(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2);
    
    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type
    operator/(const Differentiable<ScalarType, Mode>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, DiffMode Mode, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType, Mode>, OtherScalar>::Type>::type
    operator/(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType, Mode>& s2);

    template<class ScalarType, DiffMode Mode>
    inline std::ostream& operator<<(std::ostream& os, const Differentiable<ScalarType, Mode>& obj) {
        return os << obj.getValue();
    }
}

namespace std {
    template<class ScalarType, Physica::Core::DiffMode Mode>
    struct numeric_limits<Physica::Core::Differentiable<ScalarType, Mode>> : public numeric_limits<ScalarType> {};
}

#include "DifferentiableImpl/DifferentiableImpl.h"
#include "DifferentiableImpl/ElementaryFunction.h"
#include "DifferentiableImpl/ProbabilityFunction.h"
