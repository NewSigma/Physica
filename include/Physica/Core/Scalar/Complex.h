/*
 * Copyright 2020-2025 Weibo He.
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
#ifdef PHYSICA_CUDA
    #include <thrust/complex.h>
#endif
#include "Physica/Core/Scalar/Real.h" // IWYU pragma: export

namespace Physica {
    using cfloat16 = Complex<float16>;
    using cfloat32 = Complex<float32>;
    using cfloat64 = Complex<float64>;

    template<Scalar T>
    class Complex<T> : public ScalarBase<Complex<T>>, public CRCoro<Complex<T>> {
        using This = Complex<T>;
        using Base = ScalarBase<This>;
        using PacketType = BestPacket<T, 2>::Type;
    public:
        using typename Base::ScalarType;
        using typename Base::MachineType;

        constexpr static bool enableSIMD = !std::is_same<T, PacketType>::value;
    private:
        using Tm = Base::MachineType;
    private:
        T re;
        T im;
    public:
        Complex() = default;
        Complex(double _Complex x);
        __host__ __device__ Complex(Tm x);
        __host__ __device__ Complex(T re_);
        __host__ __device__ Complex(T re_, T im_);
        Complex(std::complex<Tm> x);
        Complex(MKL_Complex8 x);
        Complex(MKL_Complex16 x);
    #ifdef PHYSICA_CUDA
        __host__ __device__ Complex(thrust::complex<Tm> x);
    #endif
        template<Scalar U>
        __host__ __device__ explicit((T::Prec < U::Prec) || Diffable<U>) Complex(const U& x);
        Complex(const This&) = default;
        Complex(This&&) noexcept = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This&) = default;
        This& operator=(This&&) = default;
        [[nodiscard]] __host__ __device__ bool operator==(const This& other) const;

        template<Scalar U>
        [[nodiscard]] __host__ __device__ auto operator+(const U& x) const requires(!Diffable<U>);
        template<Scalar U>
        [[nodiscard]] __host__ __device__ auto operator-(const U& x) const requires(!Diffable<U>);
        template<Scalar U>
        [[nodiscard]] __host__ __device__ auto operator*(const U& x) const requires(!Diffable<U>);
        template<Scalar U>
        [[nodiscard]] __host__ __device__ auto operator/(const U& x) const requires(!Diffable<U>);
        [[nodiscard]] __host__ __device__ This operator-() const;
        /* Operations */
        [[nodiscard]] __host__ __device__ T squaredNorm() const;
        [[nodiscard]] __host__ __device__ T norm() const;
        [[nodiscard]] __host__ __device__ T phase() const;

        [[nodiscard]] PacketType packet() const;
        void writePacket(const PacketType packet);
        __host__ __device__ void swap(Complex& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ T& real() noexcept { return re; }
        [[nodiscard]] __host__ __device__ const T& real() const noexcept { return re; }
        [[nodiscard]] __host__ __device__ T& imag() noexcept { return im; }
        [[nodiscard]] __host__ __device__ const T& imag() const noexcept { return im; }
        [[nodiscard]] __host__ __device__ std::complex<Tm> toMachine() const noexcept;
        [[nodiscard]] __host__ __device__ auto toMachineThrust() const noexcept;
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return re.isZero() && im.isZero(); }
        [[nodiscard]] __host__ __device__ bool isFinite() const noexcept { return re.isFinite() && im.isFinite(); }
        /* Static Members */
        [[nodiscard]] static This nan() noexcept;
        [[nodiscard]] static This fromPhase(T phase) noexcept;
        template<RNG R>
        [[nodiscard]] static This random_uniform();
        template<RNG R>
        [[nodiscard]] static This random_normal();
        [[nodiscard]] static const H5::DataType& dtype_hdf5() noexcept;
    };

    template<Scalar T>
    std::ostream& operator<<(std::ostream& os, const Complex<T>& x) {
        return os << std::format("{}", x);
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<Complex<T>> {
        static_assert(!T::isComplex, "[Error]: Double complex mark is not allowed");
        static_assert(!T::isDiffable, "[Error]: Diff mark should locate in outsite");
    public:
        constexpr static FloatPrec Prec = Traits<T>::Prec;
        constexpr static int Order = 0;
        constexpr static bool isComplex = true;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;

        using ValueType = Complex<T>;
        using ScalarType = ValueType;
        using PtrTy = ScalarType*;
        using ConstPtrTy = const ScalarType*;
        using RefTy = ScalarType&;
        using ConstRefTy = const ScalarType&;
        using RealType = T;
        using ComplexType = ScalarType;
        using MachineType = T::MachineType;
    };
}

namespace std {
    template<Physica::Scalar T>
    struct numeric_limits<Physica::Complex<T>> : public numeric_limits<T> {};

    template<Physica::Scalar T>
    void swap(Physica::Complex<T>& __restrict c1, Physica::Complex<T>& __restrict c2) noexcept { c1.swap(c2); }

    template<Physica::Scalar T>
    struct formatter<Physica::Complex<T>, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        auto format(const Physica::Complex<T>& obj, std::format_context& ctx) const {
            const auto& im = obj.imag();
            return std::format_to(ctx.out(), "{} {} {}i", obj.real(), im.isNegative() ? '-' : '+', abs(im));
        }
    };
}

#include "ComplexImpl/ComplexImpl.h"
#include "ComplexImpl/SIMD.h"
