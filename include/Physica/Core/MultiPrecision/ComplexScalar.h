/*
 * Copyright 2020-2023 WeiBo He.
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

#include "Physica/Core/MultiPrecision/Scalar.h"

namespace Physica::Core {
    template<class T> class ComplexScalar;

    namespace Internal {
        template<class T>
        class Traits<ComplexScalar<T>> {
        public:
            using ScalarType = ComplexScalar<T>;
            using RealType = T;
            using ComplexType = ComplexScalar<T>;
            using TrivialType = typename T::TrivialType;
            static constexpr ScalarOption option = Traits<T>::option;
            static constexpr bool errorTrack = Traits<T>::errorTrack;
            static constexpr bool isComplex = true;
            static constexpr bool isDifferentiable = false;
        };
    }

    template<class T>
    class ComplexScalar : public ScalarBase<ComplexScalar<T>> {
        static_assert(!T::isComplex, "[Error]: Double complex mark is not allowed");
        static_assert(!T::isDifferentiable, "[Error]: Differentiable mark should locate in outsite");
        using Base = ScalarBase<ComplexScalar<T>>;
        using PacketType = typename Internal::BestPacket<T, 2>::Type;
    public:
        using typename Base::ScalarType;
        using typename Base::TrivialType;
        constexpr static bool enableSIMD = !std::is_same<T, PacketType>::value;
    private:
        T real;
        T imag;
    public:
        ComplexScalar() = default;
        ComplexScalar(T real_);
        ComplexScalar(T real_, T imag_);
        ComplexScalar(std::initializer_list<T> list);
        ComplexScalar(const ComplexScalar& c) = default;
        ComplexScalar(ComplexScalar&& c) noexcept = default;
        /* Operators */
        ComplexScalar& operator=(ComplexScalar c);
        void operator<<=(int i) { real <<= i; imag<<= i; }
        void operator>>=(int i) { real >>= i; imag>>= i; }
        bool operator==(const ComplexScalar<T>& c) const;
        bool operator!=(const ComplexScalar<T>& c) const { return !(operator==(c)); }
        /* Operations */
        [[nodiscard]] inline PacketType packet() const;
        inline void writePacket(const PacketType packet);
        void swap(ComplexScalar& c) noexcept;
        /* Getters */
        [[nodiscard]] T squaredNorm() const;
        [[nodiscard]] inline T norm() const;
        [[nodiscard]] T arg() const;
        [[nodiscard]] ComplexScalar unit() const { return *this * reciprocal(norm()); }
        [[nodiscard]] ComplexScalar conjugate() const noexcept { return ComplexScalar(real, -imag); }
        [[nodiscard]] const T& getReal() const { return real; }
        [[nodiscard]] const T& getImag() const { return imag; }
        [[nodiscard]] bool isZero() const { return real.isZero() && imag.isZero(); }
        /* Setters */
        void setReal(const T& s) { real = s; }
        void setImag(const T& s) { imag = s; }
        /* Static Members */
        template<class RandomGenerator>
        [[nodiscard]] static ComplexScalar random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static ComplexScalar random_normal(RandomGenerator& gen);
    };

    template<class T>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<T>& c);

    template<class T>
    inline ComplexScalar<T> operator+(const ComplexScalar<T>& c1, const ComplexScalar<T>& c2);

    template<class T>
    inline ComplexScalar<T> operator-(const ComplexScalar<T>& c1, const ComplexScalar<T>& c2);

    template<class T>
    ComplexScalar<T> operator*(const ComplexScalar<T>& c1, const ComplexScalar<T>& c2);

    template<class T>
    ComplexScalar<T> operator/(const ComplexScalar<T>& c1, const ComplexScalar<T>& c2);

    template<class ScalarType, ScalarOption option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<option>>::Type> operator+(
            const ComplexScalar<ScalarType>& c,const Scalar<option>& s);

    template<class ScalarType, ScalarOption option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<option>>::Type> operator-(
            const ComplexScalar<ScalarType>& c, const Scalar<option>& s);

    template<class ScalarType, ScalarOption option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<option>>::Type> operator*(
            const ComplexScalar<ScalarType>& c, const Scalar<option>& s);

    template<class ScalarType, ScalarOption option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<option>>::Type> operator/(
            const ComplexScalar<ScalarType>& c, const Scalar<option>& s);

    template<class ScalarType, ScalarOption option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<option>>::Type> operator+(
            const Scalar<option>& s, const ComplexScalar<ScalarType>& c);

    template<class ScalarType, ScalarOption option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<option>>::Type> operator-(
            const Scalar<option>& s, const ComplexScalar<ScalarType>& c);

    template<class ScalarType, ScalarOption option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<option>>::Type> operator*(
            const Scalar<option>& s, const ComplexScalar<ScalarType>& c);

    template<class ScalarType, ScalarOption option>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<option>>::Type> operator/(
            const Scalar<option>& s, const ComplexScalar<ScalarType>& c);

    template<class T>
    ComplexScalar<T> operator<<(const ComplexScalar<T>& c, int i) {
        return ComplexScalar<T>(c.getReal() << i, c.getImag() << i);
    }

    template<class T>
    ComplexScalar<T> operator>>(const ComplexScalar<T>& c, int i) {
        return ComplexScalar<T>(c.getReal() >> i, c.getImag() >> i);
    }

    template<class T>
    ComplexScalar<T> operator-(const ComplexScalar<T>& c) { return ComplexScalar(-c.getReal(), -c.getImag()); }

    template<class T>
    void swap(ComplexScalar<T>& c1, ComplexScalar<T>& c2) noexcept { c1.swap(c2); }

    template<class T, class ScalarType>
    void operator+=(ComplexScalar<T>& c, const ScalarType& t) { c = c + t; }

    template<class T, class ScalarType>
    void operator-=(ComplexScalar<T>& c, const ScalarType& t) { c = c - t; }

    template<class T, class ScalarType>
    void operator*=(ComplexScalar<T>& c, const ScalarType& t) { c = c * t; }

    template<class T, class ScalarType>
    void operator/=(ComplexScalar<T>& c, const ScalarType& t) { c = c / t; }
}

namespace std {
    template<class RealType>
    struct numeric_limits<Physica::Core::ComplexScalar<RealType>> : public numeric_limits<RealType> {};
}

#include "ComplexScalarImpl/ComplexScalarImpl.h"
