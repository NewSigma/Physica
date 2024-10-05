/*
 * Copyright 2020-2023 Weibo He.
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
    Complex<T>::Complex(double d) : Complex(T(d)) {}

    template<class T>
    Complex<T>::Complex(T real_)
            : real(real_), imag(0) {}

    template<class T>
    Complex<T>::Complex(T real_, T imag_)
            : real(real_), imag(imag_) {}

    template<class T>
    Complex<T>::Complex(std::initializer_list<T> list) {
        assert(list.size() == 2);
        auto ite = list.begin();
        real = *ite;
        ++ite;
        imag = *ite;
    }

    template<class T>
    Complex<T>& Complex<T>::operator=(Complex c) {
        swap(c);
        return *this;
    }

    template<class T>
    bool Complex<T>::operator==(const Complex<T>& c) const {
        return real == c.real && imag == c.imag;
    }

    template<class T>
    inline typename Complex<T>::PacketType Complex<T>::packet() const {
        PacketType packet{};
        packet.load(&real);
        return packet;
    }

    template<class T>
    inline void Complex<T>::writePacket(const PacketType packet) {
        packet.store(&real);
    }

    template<class T>
    void Complex<T>::swap(Complex& __restrict c) noexcept {
        assert(this != &c && "[Error]: Self swap is likely a bug");
        real.swap(c.real);
        imag.swap(c.imag);
    }

    template<class T>
    T Complex<T>::squaredNorm() const {
        return square(real) + square(imag);
    }

    template<class T>
    inline T Complex<T>::norm() const {
        return sqrt(squaredNorm());
    }

    template<class T>
    T Complex<T>::phase() const {
        auto result = arctan(imag / real);
        //Result of arctan is defined on [-Pi / 2, Pi / 2], Result of arg is defined on [-Pi, Pi].
        if(real.isNegative() && !imag.isZero())
            result += T(M_PI);
        return result;
    }

    template<class T>
    Complex<T> Complex<T>::unit() const {
        const T temp = norm();
        if (temp.isZero())
            return Complex(1);
        return *this * reciprocal(temp);
    }

    template<class T>
    inline Complex<T> Complex<T>::fromPhase(T phase) {
        T s, c;
        sincos(phase, s, c);
        return Complex(c, s);
    }

    template<class T>
    template<class RandomGenerator>
    Complex<T> Complex<T>::random_uniform(RandomGenerator& gen) {
        T real = T::random_uniform(gen);
        T imag = T::random_uniform(gen);
        return Complex(std::move(real), std::move(imag));
    }

    template<class T>
    template<class RandomGenerator>
    Complex<T> Complex<T>::random_normal(RandomGenerator& gen) {
        T real = T::random_normal(gen);
        T imag = T::random_normal(gen);
        return Complex(std::move(real), std::move(imag));
    }
#ifdef PHYSICA_HDF5
    template<class T>
    const H5::DataType& Complex<T>::getH5DataType() {
        static const auto instance = std::unique_ptr<H5::DataType>([]() -> H5::DataType* {
            auto* result = new H5::DataType(H5T_COMPOUND, sizeof(This));
            const auto id = result->getId();
            H5Tinsert(id, "Real", HOFFSET(This, real), T::getH5DataType().getId());
            H5Tinsert(id, "Imag", HOFFSET(This, imag), T::getH5DataType().getId());
            return result;
        }());
        return *instance;
    }
#endif
    template<class T>
    std::ostream& operator<<(std::ostream& os, const Complex<T>& c) {
        const auto& imag = c.getImag();
        return os << c.getReal()
                  << (imag.isNegative() ? " - " : " + " )
                  << abs(imag) << 'i';
    }

    template<class T>
    inline Complex<T> operator+(const Complex<T>& c1, const Complex<T>& c2) {
        return Complex<T>(c1.getReal() + c2.getReal(), c1.getImag() + c2.getImag());
    }

    template<class T>
    inline Complex<T> operator-(const Complex<T>& c1, const Complex<T>& c2) {
        return Complex<T>(c1.getReal() - c2.getReal(), c1.getImag() - c2.getImag());
    }

    template<class T>
    Complex<T> operator*(const Complex<T>& c1, const Complex<T>& c2) {
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
        return Complex<T>(ac - bd
                , (real_1 + imagine_1) * (real_2 + imagine_2) - ac - bd);
    }

    template<class T>
    Complex<T> operator/(const Complex<T>& c1, const Complex<T>& c2) {
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
        return Complex<T>((ac + bd) / divisor
                , ((real_1 + imagine_1) * (real_2 - imagine_2) - ac + bd) / divisor);
    }

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator+(
            const Complex<ScalarType>& c,const Scalar<Option>& s) {
        return {c.getReal() + s, c.getImag()};
    }

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator-(
            const Complex<ScalarType>& c, const Scalar<Option>& s) {
        return {c.getReal() - s, c.getImag()};
    }

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator*(
            const Complex<ScalarType>& c, const Scalar<Option>& s) {
        return {c.getReal() * s, c.getImag() * s};
    }

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator/(
            const Complex<ScalarType>& c, const Scalar<Option>& s) {
        const auto rep = reciprocal(s);
        return {c.getReal() * rep, c.getImag() * rep};
    }

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator+(
            const Scalar<Option>& s, const Complex<ScalarType>& c) {
        return c + s;
    }

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator-(
            const Scalar<Option>& s, const Complex<ScalarType>& c) {
        return {s - c.getReal(), c.getImag()};
    }

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator*(
            const Scalar<Option>& s, const Complex<ScalarType>& c) {
        return c * s;
    }

    template<class ScalarType, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpReturnType<ScalarType, Scalar<Option>>::Type> operator/(
            const Scalar<Option>& s, const Complex<ScalarType>& c) {
        const auto& real = c.getReal();
        const auto& imag = c.getImag();
        const auto divisor = s * reciprocal(square(real) + square(imag));
        return {real * divisor, -imag * divisor};
    }
}

#include "ElemFunc.h"
