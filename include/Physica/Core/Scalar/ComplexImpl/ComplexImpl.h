/*
 * Copyright 2020-2026 Weibo He.
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

#include "../Complex.h"

namespace Physica {
    template<Scalar T>
    Complex<T>::Complex(double _Complex x) : This(std::complex<Tm>(x)) {}

    template<Scalar T>
    __host__ __device__ Complex<T>::Complex(Tm x) : This(T(x)) {}

    template<Scalar T>
    __host__ __device__ Complex<T>::Complex(T re_) : re(re_), im(0) {}

    template<Scalar T>
    __host__ __device__ Complex<T>::Complex(T re_, T im_) : re(re_), im(im_) {}

    template<Scalar T>
    Complex<T>::Complex(std::complex<Tm> x) : re(x.real()), im(x.imag()) {}

#ifdef PHYSICA_MKL
    template<Scalar T>
    Complex<T>::Complex(MKL_Complex8 x) : This(x.real, x.imag) {}

    template<Scalar T>
    Complex<T>::Complex(MKL_Complex16 x) : This(x.real, x.imag) {}
#endif

#ifdef PHYSICA_CUDA
    template<Scalar T>
    __host__ __device__ Complex<T>::Complex(thrust::complex<Tm> x) : re(x.real()), im(x.imag()) {}
#endif

    template<Scalar T>
    template<Scalar U>
    __host__ __device__ Complex<T>::Complex(const U& x) {
        if constexpr (Diffable<U>)
            *this = x.value();
        else if constexpr (U::isComplex)
            *this = This(T(x.real()), T(x.imag()));
        else
            *this = This(T(x.real()));
    }

    template<Scalar T>
    Complex<T>::operator MKL_Complex() const noexcept {
        return toMKL();
    }

    template<Scalar T>
    Complex<T>::operator cuBLAS_Complex() const noexcept {
        return toCUDA();
    }

    template<Scalar T>
    __host__ __device__ bool Complex<T>::operator==(const This& other) const {
        return re == other.re && im == other.im;
    }

    template<Scalar T>
    template<Scalar U>
    __host__ __device__ auto Complex<T>::operator+(const U& x) const noexcept requires(!Diffable<U>) {
        using RtnType = Internal::BinaryScalarOpRtnTy<This, U>::Type;
        if constexpr (U::isComplex)
            return RtnType(real() + x.real(), imag() + x.imag());
        else
            return RtnType(real() + x.real(), imag());
    }

    template<Scalar T>
    template<Scalar U>
    __host__ __device__ auto Complex<T>::operator-(const U& x) const noexcept requires(!Diffable<U>) {
        using RtnType = Internal::BinaryScalarOpRtnTy<This, U>::Type;
        if constexpr (U::isComplex)
            return RtnType(real() - x.real(), imag() - x.imag());
        else
            return RtnType(real() - x.real(), imag());
    }

    template<Scalar T>
    template<Scalar U>
    __host__ __device__ auto Complex<T>::operator*(const U& x) const noexcept requires(!Diffable<U>) {
        using RtnType = Internal::BinaryScalarOpRtnTy<This, U>::Type;
        if constexpr (U::isComplex) {
            const auto& re_1 = real();
            const auto& im_1 = imag();
            const auto& re_2 = x.real();
            const auto& im_2 = x.imag();
            /*
             * Optimize:
             * Use (a + ib)(c + id) = (ac - bd) + i((a + b)(c + d) - ac - bd)
             * instead of (a + ib)(c + id) = (ac - bd) + i(ad + bc) to avoid multiply.
             * But it is unclear if this method is useful to every machine.
             * May be add checks and use Config.h to determine which method to use.
             */
            const auto ac = re_1 * re_2;
            const auto bd = im_1 * im_2;
            return RtnType(ac - bd, (re_1 + im_1) * (re_2 + im_2) - ac - bd);
        }
        else
            return RtnType(real() * x, imag() * x);
    }

    template<Scalar T>
    template<Scalar U>
    __host__ __device__ auto Complex<T>::operator/(const U& x) const noexcept requires(!Diffable<U>) {
        assert(!x.isSubNormal() && "[Error]: Division overflow");
        using RtnType = Internal::BinaryScalarOpRtnTy<This, U>::Type;
        if constexpr (U::isComplex) {
            const auto& re_1 = real();
            const auto& im_1 = imag();
            const auto& re_2 = x.real();
            const auto& im_2 = x.imag();
            /*
             * Optimize: Using the same method with operator*().
             */
            const auto ac = re_1 * re_2;
            const auto bd = im_1 * im_2;
            // May overflow
            // Algorithm 116; https://dl.acm.org/doi/pdf/10.1145/368637.368661
            const auto divisor = square(re_2) + square(im_2);
            return RtnType((ac + bd) / divisor, ((re_1 + im_1) * (re_2 - im_2) - ac + bd) / divisor);
        }
        else {
            const auto rep = reciprocal(x);
            return RtnType(real() * rep, imag() * rep);
        }
    }

    template<Scalar T>
    __host__ __device__ auto Complex<T>::operator-() const -> This {
        return Complex<T>(-real(), -imag());
    }

    template<Scalar T>
    __host__ __device__ T Complex<T>::squaredNorm() const {
        return square(re) + square(im);
    }

    template<Scalar T>
    __host__ __device__ T Complex<T>::norm() const {
        return sqrt(squaredNorm());
    }

    template<Scalar T>
    __host__ __device__ T Complex<T>::phase() const {
        return std::arg(toMachine());
    }

    template<Scalar T>
    auto Complex<T>::packet() const -> PacketType {
        PacketType packet{};
        packet.load(&re);
        return packet;
    }

    template<Scalar T>
    void Complex<T>::writePacket(const PacketType packet) {
        packet.store(&re);
    }

    template<Scalar T>
    __host__ __device__ void Complex<T>::swap(Complex& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        re.swap(obj.re);
        im.swap(obj.im);
    }

    template<Scalar T>
    __host__ __device__ auto Complex<T>::toMachine() const noexcept -> std::complex<Tm> {
        return {re.toMachine(), im.toMachine()};
    }

    template<Scalar T>
    __host__ __device__ auto Complex<T>::toThrust() const noexcept {
    #ifdef PHYSICA_CUDA
        return thrust::complex<Tm>(toMachine());
    #else
        return;
    #endif
    }

    template<Scalar T>
    __host__ __device__ auto Complex<T>::toMKL() const noexcept -> MKL_Complex {
        return {re.toMachine(), im.toMachine()};
    }

    template<Scalar T>
    __host__ __device__ auto Complex<T>::toCUDA() const noexcept -> cuBLAS_Complex {
        return {re.toMachine(), im.toMachine()};
    }

    template<Scalar T>
    auto Complex<T>::nan() noexcept -> Complex<T> {
        return {T::nan(), T::nan()};
    }

    template<Scalar T>
    auto Complex<T>::fromPhase(T phase) noexcept -> Complex<T> {
        This result{};
        sincos(phase, result.im, result.re);
        return result;
    }

    template<Scalar T>
    template<RNG R>
    auto Complex<T>::random_uniform() -> Complex<T> {
        T re = T::template random_uniform<R>();
        T im = T::template random_uniform<R>();
        return Complex(std::move(re), std::move(im));
    }

    template<Scalar T>
    template<RNG R>
    auto Complex<T>::random_normal() -> Complex<T> {
        T re = T::template random_normal<R>();
        T im = T::template random_normal<R>();
        return Complex(std::move(re), std::move(im));
    }

#ifdef PHYSICA_HDF5
    template<Scalar T>
    const H5::DataType& Complex<T>::dtype_hdf5() noexcept {
        static const auto instance = std::unique_ptr<H5::DataType>([]() -> H5::DataType* {
            auto* result = new (std::nothrow) H5::DataType(H5T_COMPOUND, sizeof(This));
            const auto id = result->getId();
            H5Tinsert(id, "Real", HOFFSET(This, re), T::dtype_hdf5().getId());
            H5Tinsert(id, "Imag", HOFFSET(This, im), T::dtype_hdf5().getId());
            return result;
        }());
        return *instance;
    }
#endif

    template<Scalar T, Scalar U>
    __host__ __device__ auto operator+(const T& x, const Complex<U>& y) requires(!T::isComplex && !Diffable<T>) {
        return y + x;
    }

    template<Scalar T, Scalar U>
    __host__ __device__ auto operator-(const T& x, const Complex<U>& y) requires(!T::isComplex && !Diffable<T>) {
        using RtnType = Internal::BinaryScalarOpRtnTy<T, Complex<U>>::Type;
        return RtnType(x - y.real(), -y.imag());
    }

    template<Scalar T, Scalar U>
    __host__ __device__ auto operator*(const T& x, const Complex<U>& y) requires(!T::isComplex && !Diffable<T>) {
        return y * x;
    }

    template<Scalar T, Scalar U>
    __host__ __device__ auto operator/(const T& x, const Complex<U>& y) requires(!T::isComplex && !Diffable<T>) {
        using RtnType = Internal::BinaryScalarOpRtnTy<T, Complex<U>>::Type;
        const auto& re = y.real();
        const auto& im = y.imag();
        const auto divisor = x * reciprocal(square(re) + square(im));
        return RtnType(re * divisor, -im * divisor);
    }
}

#include "Math.h"
