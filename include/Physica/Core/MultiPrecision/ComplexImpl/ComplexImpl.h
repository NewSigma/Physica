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

namespace Physica::Core {
    template<Scalar T>
    Complex<T>::Complex(T re_) : re(re_), im(0) {}

    template<Scalar T>
    Complex<T>::Complex(T re_, T im_) : re(re_), im(im_) {}

    template<Scalar T>
    Complex<T>::Complex(std::initializer_list<T> list) {
        assert(list.size() == 2);
        auto ite = list.begin();
        re = *ite;
        ++ite;
        im = *ite;
    }

    template<Scalar T>
    Complex<T>::Complex(std::complex<MachineType> c) : re(c.real()), im(c.imag()) {}

    template<Scalar T>
    template<Scalar U, DiffMode Mode, int Order>
    Complex<T>::Complex(const Diff<U, Mode, Order>& d) : Complex(d.value()) {}

    template<Scalar T>
    inline bool Complex<T>::operator==(const This& other) const {
        return re == other.re && im == other.im;
    }

    template<Scalar T>
    inline T Complex<T>::squaredNorm() const {
        return square(re) + square(im);
    }

    template<Scalar T>
    inline T Complex<T>::norm() const {
        return sqrt(squaredNorm());
    }

    template<Scalar T>
    inline T Complex<T>::phase() const {
        return std::arg(toMachine());
    }

    template<Scalar T>
    Complex<T> Complex<T>::unit() const {
        const T temp = norm();
        if (temp.isZero())
            return T(1);
        return *this * reciprocal(temp);
    }

    template<Scalar T>
    inline Complex<T>::PacketType Complex<T>::packet() const {
        PacketType packet{};
        packet.load(&re);
        return packet;
    }

    template<Scalar T>
    inline void Complex<T>::writePacket(const PacketType packet) {
        packet.store(&re);
    }

    template<Scalar T>
    void Complex<T>::swap(Complex& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        re.swap(obj.re);
        im.swap(obj.im);
    }

    template<Scalar T>
    inline std::complex<typename Complex<T>::MachineType> Complex<T>::toMachine() const noexcept {
        return {re.toMachine(), im.toMachine()};
    }

    template<Scalar T>
    inline Complex<T> Complex<T>::fromPhase(T phase) {
        T s, c;
        sincos(phase, s, c);
        return Complex(c, s);
    }

    template<Scalar T>
    template<RandomGenerator R>
    Complex<T> Complex<T>::random_uniform() {
        T re = T::template random_uniform<R>();
        T im = T::template random_uniform<R>();
        return Complex(std::move(re), std::move(im));
    }

    template<Scalar T>
    template<RandomGenerator R>
    Complex<T> Complex<T>::random_normal() {
        T re = T::template random_normal<R>();
        T im = T::template random_normal<R>();
        return Complex(std::move(re), std::move(im));
    }
#ifdef PHYSICA_HDF5
    template<Scalar T>
    const H5::DataType& Complex<T>::getH5DataType() {
        static const auto instance = std::unique_ptr<H5::DataType>([]() -> H5::DataType* {
            auto* result = new H5::DataType(H5T_COMPOUND, sizeof(This));
            const auto id = result->getId();
            H5Tinsert(id, "Real", HOFFSET(This, re), T::getH5DataType().getId());
            H5Tinsert(id, "Imag", HOFFSET(This, im), T::getH5DataType().getId());
            return result;
        }());
        return *instance;
    }
#endif
    template<Scalar T>
    std::ostream& operator<<(std::ostream& os, const Complex<T>& c) {
        const auto& im = c.imag();
        return os << c.real()
                  << (im.isNegative() ? " - " : " + " )
                  << abs(im) << 'i';
    }

    template<Scalar T>
    inline Complex<T> operator+(const Complex<T>& c1, const Complex<T>& c2) {
        return Complex<T>(c1.real() + c2.real(), c1.imag() + c2.imag());
    }

    template<Scalar T>
    inline Complex<T> operator-(const Complex<T>& c1, const Complex<T>& c2) {
        return Complex<T>(c1.real() - c2.real(), c1.imag() - c2.imag());
    }

    template<Scalar T>
    Complex<T> operator*(const Complex<T>& c1, const Complex<T>& c2) {
        const auto& re_1 = c1.real();
        const auto& im_1 = c1.imag();
        const auto& re_2 = c2.real();
        const auto& im_2 = c2.imag();
        /*
         * Optimize:
         * Use (a + ib)(c + id) = (ac - bd) + i((a + b)(c + d) - ac - bd)
         * instead of (a + ib)(c + id) = (ac - bd) + i(ad + bc) to avoid multiply.
         * But it is unclear if this method is useful to every machine.
         * May be add checks and use Config.h to determine which method to use.
         */
        const auto ac = re_1 * re_2;
        const auto bd = im_1 * im_2;
        return Complex<T>(ac - bd
                , (re_1 + im_1) * (re_2 + im_2) - ac - bd);
    }

    template<Scalar T>
    Complex<T> operator/(const Complex<T>& c1, const Complex<T>& c2) {
        const auto& re_1 = c1.real();
        const auto& im_1 = c1.imag();
        const auto& re_2 = c2.real();
        const auto& im_2 = c2.imag();
        /*
         * Optimize: Using the same method with operator*().
         */
        const auto ac = re_1 * re_2;
        const auto bd = im_1 * im_2;
        const auto divisor = square(re_2) + square(im_2);
        return Complex<T>((ac + bd) / divisor
                , ((re_1 + im_1) * (re_2 - im_2) - ac + bd) / divisor);
    }

    template<Scalar T, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpRtnTy<T, Real<Option>>::Type> operator+(
            const Complex<T>& c,const Real<Option>& s) {
        return {c.real() + s, c.imag()};
    }

    template<Scalar T, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpRtnTy<T, Real<Option>>::Type> operator-(
            const Complex<T>& c, const Real<Option>& s) {
        return {c.real() - s, c.imag()};
    }

    template<Scalar T, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpRtnTy<T, Real<Option>>::Type> operator*(
            const Complex<T>& c, const Real<Option>& s) {
        return {c.real() * s, c.imag() * s};
    }

    template<Scalar T, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpRtnTy<T, Real<Option>>::Type> operator/(
            const Complex<T>& c, const Real<Option>& s) {
        const auto rep = reciprocal(s);
        return {c.real() * rep, c.imag() * rep};
    }

    template<Scalar T, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpRtnTy<T, Real<Option>>::Type> operator+(
            const Real<Option>& s, const Complex<T>& c) {
        return c + s;
    }

    template<Scalar T, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpRtnTy<T, Real<Option>>::Type> operator-(
            const Real<Option>& s, const Complex<T>& c) {
        return {s - c.real(), c.imag()};
    }

    template<Scalar T, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpRtnTy<T, Real<Option>>::Type> operator*(
            const Real<Option>& s, const Complex<T>& c) {
        return c * s;
    }

    template<Scalar T, ScalarOption Option>
    Complex<typename Internal::BinaryScalarOpRtnTy<T, Real<Option>>::Type> operator/(
            const Real<Option>& s, const Complex<T>& c) {
        const auto& re = c.real();
        const auto& im = c.imag();
        const auto divisor = s * reciprocal(square(re) + square(im));
        return {re * divisor, -im * divisor};
    }
}

#include "Math.h"
