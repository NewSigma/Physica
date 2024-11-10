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
    template<class T>
    class Complex : public ScalarBase<Complex<T>> {
        using This = Complex<T>;
        using Base = ScalarBase<This>;
        using PacketType = typename BestPacket<T, 2>::Type;
    public:
        using typename Base::ScalarType;
        using typename Base::MachineType;
        constexpr static bool enableSIMD = !std::is_same<T, PacketType>::value;
    #ifdef PHYSICA_MKL
        using MKL_Complex = typename std::conditional<T::Option == Float32, MKL_Complex8, MKL_Complex16>::type;
    #endif
    private:
        T re;
        T im;
    public:
        Complex() = default;
        Complex(MachineType x) : This(T(x)) {}
        Complex(T re_);
        Complex(T re_, T im_);
        Complex(std::initializer_list<T> list);
        explicit Complex(std::complex<MachineType> c);
        template<class U, DiffMode Mode, int Order>
        explicit Complex(const Diff<U, Mode, Order>& d);
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
        [[nodiscard]] Complex conjugate() const noexcept { return Complex(re, -im); }

        [[nodiscard]] inline PacketType packet() const;
        inline void writePacket(const PacketType packet);
        void swap(Complex& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T& real() noexcept { return re; }
        [[nodiscard]] const T& real() const noexcept { return re; }
        [[nodiscard]] T& imag() noexcept { return im; }
        [[nodiscard]] const T& imag() const noexcept { return im; }
        [[nodiscard]] inline std::complex<MachineType> toMachine() const noexcept;
        [[nodiscard]] bool isZero() const noexcept { return re.isZero() && im.isZero(); }
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
    Complex<T> operator-(const Complex<T>& c) { return Complex(-c.real(), -c.imag()); }

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
        static_assert(!T::isDifferentiable, "[Error]: Diff mark should locate in outsite");
    public:
        constexpr static ScalarOption Option = Traits<T>::Option;
        constexpr static int Order = 0;
        constexpr static bool isComplex = true;
        constexpr static bool isDifferentiable = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;

        using PlainScalar = Complex<T>;
        using ScalarType = PlainScalar;
        using PtrTy = ScalarType*;
        using ConstPtrTy = const ScalarType*;
        using RefTy = ScalarType&;
        using ConstRefTy = const ScalarType&;
        using RealType = T;
        using ComplexType = ScalarType;
        using MachineType = typename T::MachineType;
    };
}

namespace std {
    template<class T>
    struct numeric_limits<Physica::Core::Complex<T>> : public numeric_limits<T> {};

    template<class T>
    void swap(Physica::Core::Complex<T>& __restrict c1, Physica::Core::Complex<T>& __restrict c2) noexcept { c1.swap(c2); }
}

#include "ComplexImpl/ComplexImpl.h"
#include "ComplexImpl/SIMD.h"
