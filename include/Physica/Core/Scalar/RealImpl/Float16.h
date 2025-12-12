/*
 * Copyright 2024-2025 Weibo He.
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

#include <format>
#include <cuda_fp16.h>
#include "../Real.h"

namespace Physica {
    template<>
    class Traits<float16> {
    public:
        constexpr static FloatPrec Prec = Float16;
        constexpr static int Order = 0;
        constexpr static bool isComplex = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;

        using ValueType = float16;
        using ScalarType = ValueType;
        using PtrTy = ScalarType*;
        using ConstPtrTy = const ScalarType*;
        using RefTy = ScalarType&;
        using ConstRefTy = const ScalarType&;
        using RealType = ScalarType;
        using ComplexType = Complex<ScalarType>;
        using MachineType = half;
    };
}

namespace Physica {
    template<>
    class Real<Float16> : public ScalarBase<Real<Float16>>, public CRCoro<Real<Float16>> {
        using This = Real<Float16>;
        using Base = ScalarBase<This>;
    public:
        using device_obj_type = This;
    private:
        half h;
    public:
        constexpr Real() = default;
        __host__ __device__ Real(half f_) noexcept : h(f_) {}
        __host__ __device__ Real(std::floating_point auto f) : h(f) {}
        __host__ __device__ Real(std::integral auto i) : h(i) {}
        template<Scalar T>
        __host__ __device__ explicit(Float16 < T::Prec) Real(const T& x) requires(!T::isComplex && !Diffable<T>);
        constexpr Real(const This&) = default;
        constexpr Real(This&&) noexcept = default;
        constexpr ~Real() = default;
        /* Operators */
        using Base::operator=;
        using Base::operator>;
        using Base::operator<;
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
        __host__ __device__ explicit operator float() const noexcept { return h; }
        __host__ __device__ explicit operator double() const noexcept { return h; }
        __host__ __device__ Real operator+(const Real& s) const noexcept { return Real(h + s.h); }
        __host__ __device__ Real operator-(const Real& s) const noexcept { return Real(h - s.h); }
        __host__ __device__ Real operator*(const Real& s) const noexcept { return Real(h * s.h); }
        __host__ __device__ Real operator/(const Real& s) const noexcept;
        __host__ __device__ Real operator-() const noexcept { return Real(-h); }
        __host__ __device__ bool operator>(const Real& s) const noexcept { return h > s.h; }
        __host__ __device__ bool operator<(const Real& s) const noexcept { return h < s.h; }
        __host__ __device__ bool operator==(const Real& s) const noexcept { return h == s.h; }
        /* Operations */
        using Base::random_uniform;
        using Base::random_normal;
        __host__ __device__ void swap(Real& __restrict s) noexcept { std::swap(h, s.h); }
        /* Getters */
        [[nodiscard]] __host__ __device__ half toMachine() const noexcept { return h; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return h == half(0); }
        [[nodiscard]] __host__ __device__ bool isPositive() const noexcept { return h > half(0); }
        [[nodiscard]] __host__ __device__ bool isNegative() const noexcept { return h < half(0); }
        [[nodiscard]] __host__ __device__ bool isFinite() const noexcept { return __hisinf(h) == 0; }
        /* Static members */
        template<RNG R>
        [[nodiscard]] static Real random_uniform() noexcept;
        template<RNG R>
        [[nodiscard]] static Real random_normal() noexcept;
        template<RNG R>
        [[nodiscard]] static Real random_any(auto& distribution) noexcept { return This(float32::random_any<R>(distribution)); }
    #ifdef PHYSICA_HDF5
        [[nodiscard]] static const H5::DataType& dtype_hdf5() noexcept { return H5::PredType::NATIVE_INT16; }
    #endif
    };

    template<Scalar T>
    __host__ __device__ Real<Float16>::Real(const T& x) requires(!T::isComplex && !Diffable<T>) : h(x.toMachine()) {}

    __host__ __device__ inline Real<Float16> Real<Float16>::operator/(const This& s) const noexcept {
        assert(!s.isZero() && "[Error]: Divide by zero");
        return Real(h / s.h);
    }

    template<RNG R>
    auto Real<Float16>::random_uniform() noexcept -> This {
        return This(float32::random_uniform<R>());
    }

    template<RNG R>
    auto Real<Float16>::random_normal() noexcept -> This {
        return This(float32::random_normal<R>());
    }

    inline std::ostream& operator<<(std::ostream& os, const Real<Float16>& s) {
        int lastPrec = static_cast<int>(os.precision());
        return os << std::setprecision(4) << float(s) << std::setprecision(lastPrec);
    }
    // Marked as __device__ because nvcc warning 20014-D, host does not have native float16 support
    [[nodiscard]] __device__ inline float16 operator ""_HF(long double x) {
        return float16(float(x));
    }

    [[nodiscard]] __device__ inline float16 operator ""_HF(unsigned long long int x) {
        return float16(int(x));
    }
}

namespace std {
    template<>
    struct PHYSICA_API formatter<Physica::Real<Physica::Float16>, char> {
        constexpr static auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
        static auto format(Physica::Real<Physica::Float16> obj, std::format_context& ctx) -> std::format_context::iterator;
    };
}
