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

namespace Physica::Core {
    namespace Internal {
        template<class T>
        class Traits<Differentiable<T>> {
            using RealT = typename T::RealType;
            using ComplexT = typename T::ComplexType;
        public:
            using ScalarType = Differentiable<T>;
            using RealType = Differentiable<RealT>;
            using ComplexType = Differentiable<ComplexT>;
            using TrivialType = typename T::TrivialType;
            static constexpr ScalarOption option = T::option;
            static constexpr bool isComplex = T::isComplex;
            static constexpr bool isDifferentiable = true;
            /* SIMD */
            using BoolSIMDType = BoolSIMD<ScalarType, 1>;
        };
    }
    /**
     * Auto differential support for scalars
     */
    template<class ScalarType>
    class Differentiable : public ScalarBase<Differentiable<ScalarType>> {
        using This = Differentiable<ScalarType>;

        ScalarType value;
        ScalarType tangent;
    public:
        Differentiable() = default;
        Differentiable(ScalarType value_);
        Differentiable(double d) : This(ScalarType(d)) {}
        Differentiable(ScalarType value_, ScalarType tangent_);
        template<class AnyScalar>
        Differentiable(const Differentiable<AnyScalar>& obj);
        Differentiable(const Differentiable&) = default;
        Differentiable(Differentiable&&) noexcept = default;
        ~Differentiable() = default;
        /* Operators */
        Differentiable& operator=(Differentiable obj) noexcept;
        [[nodiscard]] explicit operator float() const { return float(value); }
        [[nodiscard]] explicit operator double() const { return double(value); }
        [[nodiscard]] inline bool operator==(const This& other);
        [[nodiscard]] bool operator!=(const This& other) { return !this->operator==(other); }
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
    };

    template<class ScalarType>
    Differentiable<ScalarType>::Differentiable(ScalarType value_)
        : value(std::move(value_)), tangent(0) {}

    template<class ScalarType>
    Differentiable<ScalarType>::Differentiable(ScalarType value_, ScalarType tangent_)
        : value(std::move(value_)), tangent(std::move(tangent_)) {}

    template<class ScalarType>
    template<class AnyScalar>
    Differentiable<ScalarType>::Differentiable(const Differentiable<AnyScalar>& obj)
            : value(obj.getValue()), tangent(obj.getTangent()) {}

    template<class ScalarType>
    Differentiable<ScalarType>& Differentiable<ScalarType>::operator=(Differentiable obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType>
    inline bool Differentiable<ScalarType>::operator==(const This& other) {
        return value == other.value && tangent == other.tangent;
    }

    template<class ScalarType>
    inline Differentiable<ScalarType> Differentiable<ScalarType>::operator-() const {
        return {-value, -tangent};
    }

    template<class ScalarType>
    void Differentiable<ScalarType>::swap(Differentiable& obj) noexcept {
        value.swap(obj.value);
        tangent.swap(obj.tangent);
    }

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType>, OtherScalar>::Type
    operator+(const Differentiable<ScalarType>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType>, OtherScalar>::Type>::type
    operator+(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType>& s2);

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType>, OtherScalar>::Type
    operator-(const Differentiable<ScalarType>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType>, OtherScalar>::Type>::type
    operator-(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType>& s2);

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType>, OtherScalar>::Type
    operator*(const Differentiable<ScalarType>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType>, OtherScalar>::Type>::type
    operator*(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType>& s2);
    
    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType>, OtherScalar>::Type
    operator/(const Differentiable<ScalarType>& s1, const ScalarBase<OtherScalar>& s2);

    template<class ScalarType, class OtherScalar>
    [[nodiscard]] inline typename std::enable_if<!OtherScalar::isDifferentiable, typename Internal::BinaryScalarOpReturnType<Differentiable<ScalarType>, OtherScalar>::Type>::type
    operator/(const ScalarBase<OtherScalar>& s1, const Differentiable<ScalarType>& s2);

    template<class ScalarType>
    inline std::ostream& operator<<(std::ostream& os, const Differentiable<ScalarType>& obj) {
        return os << obj.getValue();
    }
}

namespace std {
    template<class ScalarType>
    struct numeric_limits<Physica::Core::Differentiable<ScalarType>> : public numeric_limits<ScalarType> {};
}

#include "DifferentiableImpl/DifferentiableImpl.h"
#include "DifferentiableImpl/ElementaryFunction.h"
#include "DifferentiableImpl/ProbabilityFunction.h"
