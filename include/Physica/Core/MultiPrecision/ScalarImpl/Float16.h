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
        using ScalarType = Core::float16;
        using RealType = ScalarType;
        using ComplexType = Core::ComplexScalar<ScalarType>;
        using TrivialType = __half;
        using PlainScalar = ScalarType;
        constexpr static Core::ScalarOption Option = Core::Float16;
        constexpr static bool isComplex = false;
        constexpr static bool isDifferentiable = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;
        constexpr static unsigned int Order = 0;
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
        __half f;
    public:
        Scalar() = default;
        __host__ __device__ Scalar(int i) : f(i) {}
        __host__ __device__ Scalar(__half f_) : f(f_) {}
        __host__ __device__ Scalar(float f_) : f(f_) {}
        __host__ __device__ Scalar(double d_) : f(d_) {}
        template<class OtherScalar>
        __host__ __device__ explicit inline Scalar(const ScalarBase<OtherScalar>& s);
        Scalar(const This&) = default;
        Scalar(This&&) noexcept = default;
        //~Scalar() = default; /* Dynamic parallelism of CUDA 12.1 does not recognize that PlainStruct is trivial */
        /* Operators */
        using Base::operator>;
        using Base::operator<;
        Scalar& operator=(const Scalar& obj) = default;
        Scalar& operator=(Scalar&& obj) noexcept = default;
        __host__ __device__ explicit operator float() const { return f; }
        __host__ __device__ explicit operator double() const { return f; }
        __host__ __device__ Scalar operator+(const Scalar& s) const { return Scalar(f + s.f); }
        __host__ __device__ Scalar operator-(const Scalar& s) const { return Scalar(f - s.f); }
        __host__ __device__ Scalar operator*(const Scalar& s) const { return Scalar(f * s.f); }
        __host__ __device__ Scalar operator/(const Scalar& s) const { return Scalar(f / s.f); }
        __host__ __device__ Scalar operator-() const noexcept { return Scalar(-f); }
        __host__ __device__ bool operator>(const Scalar& s) const { return f > s.f; }
        __host__ __device__ bool operator<(const Scalar& s) const { return f < s.f; }
        __host__ __device__ bool operator==(const Scalar& s) const { return f == s.f; }
        /* Operations */
        void swap(Scalar& __restrict s) noexcept { std::swap(f, s.f); }
        /* Getters */
        [[nodiscard]] constexpr static ScalarOption getOption() { return Float16; }
        [[nodiscard]] __host__ __device__ __half getTrivial() const noexcept { return f; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return f == __half(0); }
        [[nodiscard]] __host__ __device__ bool isPositive() const noexcept { return f > __half(0); }
        [[nodiscard]] __host__ __device__ bool isNegative() const noexcept { return f < __half(0); }
        [[nodiscard]] __host__ __device__ bool isFinite() const noexcept { return __hisinf(f) == 0; }
    };

    template<class OtherScalar>
    __host__ __device__ inline Scalar<Float16>::Scalar(const ScalarBase<OtherScalar>& s) : f(s.getDerived().getTrivial()) {}

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

#include "FunctionImpl/MathFP16.h"
