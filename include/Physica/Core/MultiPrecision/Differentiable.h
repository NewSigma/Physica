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

namespace Physica::Core {
    namespace Internal {
        template<class T, DiffMode M>
        class Traits<Differentiable<T, M>> {
            static_assert(!Traits<T>::isDifferentiable, "[Error]: Nested Differentiable is not allowed");
            using RealT = typename T::RealType;
            using ComplexT = typename T::ComplexType;
        public:
            using ScalarType = Differentiable<T, M>;
            using RealType = Differentiable<RealT, M>;
            using ComplexType = Differentiable<ComplexT, M>;
            using TrivialType = typename T::TrivialType;
            using PlainScalar = T;
            constexpr static ScalarOption option = T::option;
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
     * Auto differential support for scalars
     */
    template<class ScalarType, DiffMode Mode>
    class Differentiable : public ScalarBase<Differentiable<ScalarType, Mode>> {
        using This = Differentiable<ScalarType, Mode>;

        ScalarType value;
        ScalarType tangent;
    public:
        Differentiable() = default;
        Differentiable(ScalarType value_);
        Differentiable(double d) : This(ScalarType(d)) {}
        Differentiable(ScalarType value_, ScalarType tangent_);
        template<class AnyScalar>
        Differentiable(const Differentiable<AnyScalar, Mode>& obj);
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

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode>::Differentiable(ScalarType value_)
        : value(std::move(value_)), tangent(0) {}

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode>::Differentiable(ScalarType value_, ScalarType tangent_)
        : value(std::move(value_)), tangent(std::move(tangent_)) {}

    template<class ScalarType, DiffMode Mode>
    template<class AnyScalar>
    Differentiable<ScalarType, Mode>::Differentiable(const Differentiable<AnyScalar, Mode>& obj)
            : value(obj.getValue()), tangent(obj.getTangent()) {}

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode>& Differentiable<ScalarType, Mode>::operator=(Differentiable obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, DiffMode Mode>
    inline bool Differentiable<ScalarType, Mode>::operator==(const This& other) {
        return value == other.value && tangent == other.tangent;
    }

    template<class ScalarType, DiffMode Mode>
    inline Differentiable<ScalarType, Mode> Differentiable<ScalarType, Mode>::operator-() const {
        return {-value, -tangent};
    }

    template<class ScalarType, DiffMode Mode>
    void Differentiable<ScalarType, Mode>::swap(Differentiable& obj) noexcept {
        value.swap(obj.value);
        tangent.swap(obj.tangent);
    }

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
