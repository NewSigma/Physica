/*
 * Copyright 2020-2024 Weibo He.
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

#include <complex>
#include <Physica/Core/MultiPrecision/Scalar.h>

namespace Physica::Core {
    template<class T> class Complex;

    namespace Internal {
        template<class AnyScalar1, class AnyScalar2>
        class BinaryScalarOpReturnType<Complex<AnyScalar1>, AnyScalar2> {
            static_assert(!AnyScalar1::isDifferentiable && !AnyScalar2::isDifferentiable, "[Error]: This class applies to plain scalar only");
        public:
            using Type = Complex<typename BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>;
        };

        template<class AnyScalar1, class AnyScalar2>
        class BinaryScalarOpReturnType<AnyScalar1, Complex<AnyScalar2>> {
            static_assert(!AnyScalar1::isDifferentiable && !AnyScalar2::isDifferentiable, "[Error]: This class applies to plain scalar only");
        public:
            using Type = Complex<typename BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>;
        };

        template<class AnyScalar1, class AnyScalar2>
        class BinaryScalarOpReturnType<Complex<AnyScalar1>, Complex<AnyScalar2>> {
            static_assert(!AnyScalar1::isDifferentiable && !AnyScalar2::isDifferentiable, "[Error]: This class applies to plain scalar only");
        public:
            using Type = Complex<typename BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>;
        };
    }

    template<class T>
    class Complex : public ScalarBase<Complex<T>> {
        using This = Complex<T>;
        using Base = ScalarBase<This>;
        using PacketType = typename BestPacket<T, 2>::Type;
    public:
        using typename Base::ScalarType;
        using typename Base::TrivialType;
        constexpr static bool enableSIMD = !std::is_same<T, PacketType>::value;
    #ifdef PHYSICA_MKL
        using MKL_Complex = typename std::conditional<T::Option == Float32, MKL_Complex8, MKL_Complex16>::type;
    #endif
    private:
        T real;
        T imag;
    public:
        Complex() = default;
        Complex(double d);
        Complex(T real_);
        Complex(T real_, T imag_);
        Complex(std::initializer_list<T> list);
        explicit Complex(std::complex<TrivialType> c);
        Complex(const This&) = default;
        Complex(This&&) noexcept = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        bool operator==(const This& other) const;
        /* Operations */
        [[nodiscard]] inline T squaredNorm() const;
        [[nodiscard]] inline T norm() const;
        [[nodiscard]] inline T phase() const;
        [[nodiscard]] Complex unit() const;
        [[nodiscard]] Complex conjugate() const noexcept { return Complex(real, -imag); }

        [[nodiscard]] inline PacketType packet() const;
        inline void writePacket(const PacketType packet);
        void swap(Complex& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T& getReal() noexcept { return real; }
        [[nodiscard]] const T& getReal() const noexcept { return real; }
        [[nodiscard]] T& getImag() noexcept { return imag; }
        [[nodiscard]] const T& getImag() const noexcept { return imag; }
        [[nodiscard]] inline std::complex<TrivialType> getTrivial() const noexcept;
        [[nodiscard]] bool isZero() const noexcept { return real.isZero() && imag.isZero(); }
        /* Static Members */
        [[nodiscard]] inline static Complex fromPhase(T phase);
        template<class RandomGenerator>
        [[nodiscard]] static Complex random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static Complex random_normal(RandomGenerator& gen);
        [[nodiscard]] static const H5::DataType& getH5DataType();
    };

    template<class T>
    std::ostream& operator<<(std::ostream& os, const Complex<T>& c);

    template<class T>
    inline Complex<T> operator+(const Complex<T>& c1, const Complex<T>& c2);

    template<class T>
    inline Complex<T> operator-(const Complex<T>& c1, const Complex<T>& c2);

    template<class T>
    Complex<T> operator*(const Complex<T>& c1, const Complex<T>& c2);

    template<class T>
    Complex<T> operator/(const Complex<T>& c1, const Complex<T>& c2);

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator+(
            const Complex<ScalarType>& c,const Scalar<Option>& s);

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator-(
            const Complex<ScalarType>& c, const Scalar<Option>& s);

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator*(
            const Complex<ScalarType>& c, const Scalar<Option>& s);

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator/(
            const Complex<ScalarType>& c, const Scalar<Option>& s);

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator+(
            const Scalar<Option>& s, const Complex<ScalarType>& c);

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator-(
            const Scalar<Option>& s, const Complex<ScalarType>& c);

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator*(
            const Scalar<Option>& s, const Complex<ScalarType>& c);

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator/(
            const Scalar<Option>& s, const Complex<ScalarType>& c);

    template<class T>
    Complex<T> operator-(const Complex<T>& c) { return Complex(-c.getReal(), -c.getImag()); }

    template<class T, class ScalarType>
    void operator+=(Complex<T>& c, const ScalarType& t) { c = c + t; }

    template<class T, class ScalarType>
    void operator-=(Complex<T>& c, const ScalarType& t) { c = c - t; }

    template<class T, class ScalarType>
    void operator*=(Complex<T>& c, const ScalarType& t) { c = c * t; }

    template<class T, class ScalarType>
    void operator/=(Complex<T>& c, const ScalarType& t) { c = c / t; }
}

namespace Physica {
    using namespace Core;

    template<class T>
    class Traits<Complex<T>> {
        static_assert(!T::isComplex, "[Error]: Double complex mark is not allowed");
        static_assert(!T::isDifferentiable, "[Error]: Differentiable mark should locate in outsite");
    public:
        using ScalarType = Complex<T>;
        using RealType = T;
        using ComplexType = Complex<T>;
        using TrivialType = typename T::TrivialType;
        using PlainScalar = ScalarType;
        constexpr static ScalarOption Option = Traits<T>::Option;
        constexpr static bool errorTrack = Traits<T>::errorTrack;
        constexpr static bool isComplex = true;
        constexpr static bool isDifferentiable = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;
        constexpr static unsigned int Order = 0;
    };
}

namespace std {
    template<class T>
    struct numeric_limits<Physica::Core::Complex<T>> : public numeric_limits<T> {};

    template<class T>
    void swap(Physica::Core::Complex<T>& __restrict c1, Physica::Core::Complex<T>& __restrict c2) noexcept { c1.swap(c2); }
}

#include "ComplexImpl/ComplexImpl.h"
