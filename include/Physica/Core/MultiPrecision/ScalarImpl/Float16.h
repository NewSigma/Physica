/*
 * Copyright 2024 Weibo He.
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

#include <cuda_fp16.h>
#include <iomanip>

namespace Physica {
    template<>
    class Traits<Core::float16> {
    public:
        constexpr static Core::ScalarOption Option = Core::Float16;
        constexpr static int Order = 0;
        constexpr static bool isComplex = false;
        constexpr static bool isDifferentiable = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;

        using PlainScalar = Core::float16;
        using ScalarType = PlainScalar;
        using PtrTy = ScalarType*;
        using ConstPtrTy = const ScalarType*;
        using RefTy = ScalarType&;
        using ConstRefTy = const ScalarType&;
        using RealType = ScalarType;
        using ComplexType = Core::Complex<ScalarType>;
        using MachineType = __half;
    };
}

namespace Physica::Core {
    template<>
    class Scalar<Float16> : public ScalarBase<Scalar<Float16>> {
        using This = Scalar<Float16>;
        using Base = ScalarBase<This>;
    public:
        using device_obj_type = This;
    private:
        __half h;
    public:
        Scalar() = default;
        __host__ __device__ Scalar(int i) : h(i) {}
        __host__ __device__ Scalar(__half f_) : h(f_) {}
        __host__ __device__ Scalar(float f_) : h(f_) {}
        __host__ __device__ Scalar(double d_) : h(d_) {}
        template<class OtherScalar>
        __host__ __device__ explicit inline Scalar(const ScalarBase<OtherScalar>& s);
        Scalar(const This&) = default;
        Scalar(This&&) noexcept = default;
        //~Scalar() = default; /* Dynamic parallelism of CUDA 12.1 does not recognize that PlainStruct is trivial */
        /* Operators */
        using Base::operator>;
        using Base::operator<;
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
        __host__ __device__ explicit operator float() const { return h; }
        __host__ __device__ explicit operator double() const { return h; }
        __host__ __device__ Scalar operator+(const Scalar& s) const { return Scalar(h + s.h); }
        __host__ __device__ Scalar operator-(const Scalar& s) const { return Scalar(h - s.h); }
        __host__ __device__ Scalar operator*(const Scalar& s) const { return Scalar(h * s.h); }
        __host__ __device__ Scalar operator/(const Scalar& s) const { return Scalar(h / s.h); }
        __host__ __device__ Scalar operator-() const noexcept { return Scalar(-h); }
        __host__ __device__ bool operator>(const Scalar& s) const { return h > s.h; }
        __host__ __device__ bool operator<(const Scalar& s) const { return h < s.h; }
        __host__ __device__ bool operator==(const Scalar& s) const { return h == s.h; }
        /* Operations */
        void swap(Scalar& __restrict s) noexcept { std::swap(h, s.h); }
        /* Getters */
        [[nodiscard]] constexpr static ScalarOption getOption() { return Float16; }
        [[nodiscard]] __host__ __device__ __half toMachine() const noexcept { return h; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return h == __half(0); }
        [[nodiscard]] __host__ __device__ bool isPositive() const noexcept { return h > __half(0); }
        [[nodiscard]] __host__ __device__ bool isNegative() const noexcept { return h < __half(0); }
        [[nodiscard]] __host__ __device__ bool isFinite() const noexcept { return __hisinf(h) == 0; }
    };

    template<class OtherScalar>
    __host__ __device__ inline Scalar<Float16>::Scalar(const ScalarBase<OtherScalar>& s) : h(s.getDerived().toMachine()) {}

    inline std::ostream& operator<<(std::ostream& os, const Scalar<Float16>& s) {
        const auto lastPrec = os.precision();
        return os << std::setprecision(4) << float(s) << std::setprecision(lastPrec);
    }

    [[nodiscard]] __host__ __device__ inline float16 operator ""_HF(long double x) {
        return float16(float(x));
    }

    [[nodiscard]] __host__ __device__ inline float16 operator ""_HF(unsigned long long int x) {
        return float16(int(x));
    }
}

#include "MathFP16.h"
