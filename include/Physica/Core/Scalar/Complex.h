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

#include <complex>
#ifdef PHYSICA_CUDA
    #include <cuComplex.h>
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
    #ifndef PHYSICA_MKL
        struct MKL_Complex8 {};
        struct MKL_Complex16 {};
    #endif

    #ifndef PHYSICA_CUDA
        struct cuComplex {};
        struct cuDoubleComplex {};
    #endif

    public:
        using typename Base::ScalarType;
        using typename Base::MachineType;
        using MKL_Complex = std::conditional<T::Prec == Float32, MKL_Complex8, typename std::conditional<T::Prec == Float64, MKL_Complex16, void>::type>::type;
        using cuBLAS_Complex = std::conditional<T::Prec == Float32, cuComplex, typename std::conditional<T::Prec == Float64, cuDoubleComplex, void>::type>::type;
    private:
        using Tm = Base::MachineType;
    private:
        T re;
        T im;
    public:
        Complex() = default;
        __host__ __device__ Complex(double _Complex x);
        __host__ __device__ Complex(Tm x);
        __host__ __device__ Complex(T re_);
        __host__ __device__ Complex(T re_, T im_);
        __host__ __device__ Complex(std::complex<Tm> x);
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
        [[nodiscard]] explicit operator MKL_Complex() const noexcept;
        [[nodiscard]] explicit operator cuBLAS_Complex() const noexcept;
        [[nodiscard]] __host__ __device__ bool operator==(const This& other) const;

        template<Scalar U>
        [[nodiscard]] __host__ __device__ auto operator+(const U& x) const noexcept requires(!Diffable<U>);
        template<Scalar U>
        [[nodiscard]] __host__ __device__ auto operator-(const U& x) const noexcept requires(!Diffable<U>);
        template<Scalar U>
        [[nodiscard]] __host__ __device__ auto operator*(const U& x) const noexcept requires(!Diffable<U>);
        template<Scalar U>
        [[nodiscard]] __host__ __device__ auto operator/(const U& x) const noexcept requires(!Diffable<U>);
        [[nodiscard]] __host__ __device__ This operator-() const;
        /* Operations */
        [[nodiscard]] __host__ __device__ T squaredNorm() const;
        [[nodiscard]] __host__ __device__ T norm() const;
        [[nodiscard]] __host__ __device__ T phase() const;

        using Base::random_uniform;
        using Base::random_normal;
        __host__ __device__ void swap(Complex& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ T& real() noexcept { return re; }
        [[nodiscard]] __host__ __device__ const T& real() const noexcept { return re; }
        [[nodiscard]] __host__ __device__ T& imag() noexcept { return im; }
        [[nodiscard]] __host__ __device__ const T& imag() const noexcept { return im; }
        [[nodiscard]] __host__ __device__ auto* real_ptr(this auto&& self) noexcept { return &self.re; }
        [[nodiscard]] __host__ __device__ auto* imag_ptr(this auto&& self) noexcept { return &self.im; }
        [[nodiscard]] __host__ __device__ std::complex<Tm> toMachine() const noexcept;
        [[nodiscard]] __host__ __device__ auto toThrust() const noexcept;
        [[nodiscard]] __host__ __device__ MKL_Complex toMKL() const noexcept;
        [[nodiscard]] __host__ __device__ cuBLAS_Complex toCUDA() const noexcept;
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return re.isZero() && im.isZero(); }
        [[nodiscard]] __host__ __device__ bool isSubNormal() const noexcept { return re.isSubNormal() && im.isSubNormal(); }
        [[nodiscard]] __host__ __device__ bool isFinite() const noexcept { return re.isFinite() && im.isFinite(); }
        /* Static Members */
        [[nodiscard]] static This nan() noexcept;
        [[nodiscard]] static This fromPhase(T phase) noexcept;
        template<RNG R>
        [[nodiscard]] static This random_uniform();
        template<RNG R>
        [[nodiscard]] static This random_normal();
        [[nodiscard]] static const H5Type& dtype_hdf5() noexcept;
    };

    template<Scalar T>
    std::ostream& operator<<(std::ostream& os, const Complex<T>& x) {
        return os << std::format("{}", x);
    }

    template<Scalar T>
    void swap(Complex<T>& c1, Complex<T>& c2) noexcept { c1.swap(c2); }
}

namespace Physica {
    template<Scalar T>
    class Traits<Complex<T>> {
        static_assert(!T::isComplex(), "[Error]: Double complex mark is not allowed");
        static_assert(!T::isDiffable(), "[Error]: Diff mark should locate in outsite");
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
    struct formatter<Physica::Complex<T>, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        auto format(const Physica::Complex<T>& obj, auto& ctx) const {
            const auto& im = obj.imag();
            return std::format_to(ctx.out(), "{} {} {}i", obj.real(), im.isNegative() ? '-' : '+', abs(im));
        }
    };
}

#include "ComplexImpl/ComplexImpl.h"
#include "ComplexImpl/SIMD.h"
