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

namespace Physica::Core {
    template<class T>
    ComplexScalar<T>::ComplexScalar(const T& real_)
            : real(real_), imag(T::Zero()) {}

    template<class T>
    ComplexScalar<T>::ComplexScalar(const T& real_, const T& imag_)
            : real(real_), imag(imag_) {}

    template<class T>
    ComplexScalar<T>::ComplexScalar(std::initializer_list<T> list) {
        assert(list.size() == 2);
        auto ite = list.begin();
        real = *ite;
        ++ite;
        imag = *ite;
    }

    template<class T>
    ComplexScalar<T>& ComplexScalar<T>::operator=(ComplexScalar c) {
        swap(c);
        return *this;
    }

    template<class T>
    bool ComplexScalar<T>::operator==(const ComplexScalar<T>& c) const {
        return real == c.real && imag == c.imag;
    }

    template<class T>
    inline typename ComplexScalar<T>::PacketType ComplexScalar<T>::packet() const {
        PacketType packet{};
        packet.load(reinterpret_cast<TrivialType*>(const_cast<T*>(&real)));
        return packet;
    }

    template<class T>
    inline void ComplexScalar<T>::writePacket(const PacketType packet) {
        packet.store(reinterpret_cast<TrivialType*>(&real));
    }

    template<class T>
    void ComplexScalar<T>::swap(ComplexScalar& c) noexcept {
        real.swap(c.real);
        imag.swap(c.imag);
    }

    template<class T>
    T ComplexScalar<T>::squaredNorm() const {
        return square(real) + square(imag);
    }

    template<class T>
    inline T ComplexScalar<T>::norm() const {
        return sqrt(squaredNorm());
    }

    template<class T>
    T ComplexScalar<T>::arg() const {
        auto result = arctan(imag / real);
        //Result of arctan is defined on [-Pi / 2, Pi / 2], Result of arg is defined on [-Pi, Pi].
        if(real.isNegative() && !imag.isZero())
            result += T(M_PI);
        return result;
    }

    template<class T>
    inline ComplexScalar<T> ComplexScalar<T>::Zero() {
        return ComplexScalar(T::Zero(), T::Zero());
    }

    template<class T>
    inline ComplexScalar<T> ComplexScalar<T>::One() {
        return ComplexScalar(T::One(), T::Zero());
    }

    template<class T>
    inline ComplexScalar<T> ComplexScalar<T>::Two() {
        return ComplexScalar(T::Two(), T::Zero());
    }

    template<class T>
    template<class RandomGenerator>
    ComplexScalar<T> ComplexScalar<T>::random_uniform(RandomGenerator& gen) {
        T real = T::random_uniform(gen);
        T imag = T::random_uniform(gen);
        return ComplexScalar(std::move(real), std::move(imag));
    }

    template<class T>
    template<class RandomGenerator>
    ComplexScalar<T> ComplexScalar<T>::random_normal(RandomGenerator& gen) {
        T real = T::random_normal(gen);
        T imag = T::random_normal(gen);
        return ComplexScalar(std::move(real), std::move(imag));
    }

    template<class T>
    bool scalarNear(const ComplexScalar<T>& s1,
                    const ComplexScalar<T>& s2,
                    double precision) {
        assert(precision > 0);
        return scalarNear((s1 - s2).norm(), T(0), precision);
    }

    template<class T>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<T>& c) {
        const auto& imagine = c.getImag();
        return os << std::setprecision(10) << double(c.getReal())
                  << (imagine.isNegative() ? " - " : " + " ) << fabs(double(imagine)) << 'i' << std::setprecision(6);
    }

    template<class T>
    ComplexScalar<T> operator+(
            const ComplexScalar<T>& c1, const ComplexScalar<T>& c2) {
        return ComplexScalar<T>(c1.getReal() + c2.getReal(), c1.getImag() + c2.getImag());
    }

    template<class T>
    ComplexScalar<T> operator-(
            const ComplexScalar<T>& c1, const ComplexScalar<T>& c2) {
        return ComplexScalar<T>(c1.getReal() - c2.getReal(), c1.getImag() - c2.getImag());
    }

    template<class T>
    ComplexScalar<T> operator*(
            const ComplexScalar<T>& c1, const ComplexScalar<T>& c2) {
        const auto& real_1 = c1.getReal();
        const auto& imagine_1 = c1.getImag();
        const auto& real_2 = c2.getReal();
        const auto& imagine_2 = c2.getImag();
        /*
         * Optimize:
         * Use (a + ib)(c + id) = (ac - bd) + i((a + b)(c + d) - ac - bd)
         * instead of (a + ib)(c + id) = (ac - bd) + i(ad + bc) to avoid multiply.
         * But it is unclear if this method is useful to every machine.
         * May be add checks and use Config.h to determine which method to use.
         */
        const auto ac = real_1 * real_2;
        const auto bd = imagine_1 * imagine_2;
        return ComplexScalar<T>(ac - bd
                , (real_1 + imagine_1) * (real_2 + imagine_2) - ac - bd);
    }

    template<class T>
    ComplexScalar<T> operator/(
            const ComplexScalar<T>& c1, const ComplexScalar<T>& c2) {
        const auto& real_1 = c1.getReal();
        const auto& imagine_1 = c1.getImag();
        const auto& real_2 = c2.getReal();
        const auto& imagine_2 = c2.getImag();
        /*
         * Optimize: Using the same method with operator*().
         */
        const auto ac = real_1 * real_2;
        const auto bd = imagine_1 * imagine_2;
        const auto divisor = square(real_2) + square(imagine_2);
        return ComplexScalar<T>((ac + bd) / divisor
                , ((real_1 + imagine_1) * (real_2 - imagine_2) - ac + bd) / divisor);
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator+(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(c.getReal() + s.getDerived(), c.getImag());
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator-(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(c.getReal() - s.getDerived(), c.getImag());
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator*(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(c.getReal() * s.getDerived(), c.getImag() * s.getDerived());
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator/(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(c.getReal() / s.getDerived(), c.getImag() / s.getDerived());
    }

    template<class ScalarType1, class ScalarType2>
    inline ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator+(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c) {
        return c + s;
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator-(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(s.getDerived() - c.getReal(), c.getImag());
    }

    template<class ScalarType1, class ScalarType2>
    inline ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator*(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c) {
        return c * s;
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator/(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        const auto divisor = s.getDerived() * reciprocal(square(real) + square(imagine));
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(real * divisor, -imagine * divisor);
    }
}

#include "CElementaryFunction.h"
