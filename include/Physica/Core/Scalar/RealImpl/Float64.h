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
#include "Physica/CRCoro.h"
#include "../Real.h"

namespace Physica {
    template<>
    class Traits<float64> {
    public:
        constexpr static FloatPrec Prec = Float64;
        constexpr static int Order = 0;
        constexpr static bool isComplex = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;

        using ValueType = float64;
        using ScalarType = ValueType;
        using PtrTy = ScalarType*;
        using ConstPtrTy = const ScalarType*;
        using RefTy = ScalarType&;
        using ConstRefTy = const ScalarType&;
        using RealType = ScalarType;
        using ComplexType = Complex<ScalarType>;
        using MachineType = double;
    };
}

namespace Physica {
    template<>
    class Real<Float64> : public ScalarBase<Real<Float64>>, public CRCoro<Real<Float64>> {
        using This = Real<Float64>;
        using Base = ScalarBase<This>;
    public:
        using device_obj_type = This;
    private:
        double d;
    public:
        constexpr Real() = default;
        __host__ __device__ constexpr Real(double d_) noexcept : d(d_) {}
        Real(const Integer& i) : Real(double(i)) {}
        Real(const Rational& r) : Real(double(r)) {}
        template<Scalar T>
        __host__ __device__ explicit(Float64 < T::Prec) Real(const T& x) requires(!T::isComplex && !Diffable<T>);
        constexpr Real(const Real&) = default;
        constexpr Real(Real&&) noexcept = default;
        constexpr ~Real() = default;
        /* Operators */
        using Base::operator=;
        using Base::operator>;
        using Base::operator<;
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
        __host__ __device__ explicit operator float() const noexcept { return d; }
        __host__ __device__ explicit operator double() const noexcept { return d; }
        __host__ __device__ Real operator+(const Real& s) const noexcept { return Real(d + s.d); }
        __host__ __device__ Real operator-(const Real& s) const noexcept { return Real(d - s.d); }
        __host__ __device__ Real operator*(const Real& s) const noexcept { return Real(d * s.d); }
        __host__ __device__ Real operator/(const Real& s) const noexcept { return Real(d / s.d); }
        __host__ __device__ Real operator<<(int i) const { return Real(std::ldexp(d, i)); }
        __host__ __device__ Real operator>>(int i) const { return Real(std::ldexp(d, -i)); }
        __host__ __device__ Real operator-() const noexcept { return Real(-d); }
        __host__ __device__ bool operator>(const Real& s) const noexcept { return d > s.d; }
        __host__ __device__ bool operator<(const Real& s) const noexcept { return d < s.d; }
        __host__ __device__ bool operator==(const Real& s) const noexcept { return d == s.d; }
        PHYSICA_API friend std::istream& operator>>(std::istream& is, Real& scalar);
        /* Operations */
        [[nodiscard]] __host__ __device__ inline Real mod() const noexcept;
        __host__ __device__ void swap(Real& __restrict s) noexcept { std::swap(d, s.d); }
        /* Getters */
        [[nodiscard]] __host__ __device__ double toMachine() const noexcept { return d; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept{ return d == 0; }
        [[nodiscard]] __host__ __device__ bool isPositive() const noexcept { return d > 0; }
        [[nodiscard]] __host__ __device__ bool isNegative() const noexcept { return d < 0; }
        [[nodiscard]] __host__ __device__ inline bool isFinite() const noexcept;
        /* Static Members */
        [[nodiscard]] inline static Real nan() noexcept;
        template<RNG R>
        [[nodiscard]] inline static Real random_uniform();
        template<RNG R>
        [[nodiscard]] inline static Real random_normal();
        template<RNG R>
        [[nodiscard]] static Real random_normal(GaussRandomPool<This, R>& pool) { return pool(); }
        template<RNG R, class Distribution>
        [[nodiscard]] inline static Real random_any(Distribution& dist);
    #ifdef PHYSICA_HDF5
        [[nodiscard]] static const H5::DataType& getH5DataType() { return H5::PredType::NATIVE_DOUBLE; }
    #endif
    };

    template<Scalar T>
    __host__ __device__ Real<Float64>::Real(const T& x) requires(!T::isComplex && !Diffable<T>) : d(double(x)) {}

    __host__ __device__ inline Real<Float64> Real<Float64>::mod() const noexcept {
        double buffer;
        return std::modf(toMachine(), &buffer);
    }

    __host__ __device__ inline bool Real<Float64>::isFinite() const noexcept {
        using namespace std;
        return isfinite(d);
    }

    inline Real<Float64> Real<Float64>::nan() noexcept {
        return std::nan("");
    }

    template<RNG R>
    inline Real<Float64> Real<Float64>::random_uniform() {
        return Real(std::generate_canonical<double, std::numeric_limits<double>::digits>(R::getInstance()));
    }

    template<RNG R>
    inline Real<Float64> Real<Float64>::random_normal() {
        std::normal_distribution<double> dist{};
        return Real(dist(R::getInstance()));
    }

    template<RNG R, class Distribution>
    inline Real<Float64> Real<Float64>::random_any(Distribution& dist) {
        return Real(dist(R::getInstance()));
    }

    inline std::ostream& operator<<(std::ostream& os, const Real<Float64>& x) {
        return os << std::format("{}", x.toMachine());
    }

    inline std::istream& operator>>(std::istream& is, Real<Float64>& scalar) {
        is >> scalar.d;
        return is;
    }
}

namespace std {
    template<>
    struct numeric_limits<Physica::Real<Physica::Float64>> : 
    #ifdef PHYSICA_CUDA
        public ::cuda::std::numeric_limits<double>
    #else
        public numeric_limits<double>
    #endif
    {};

    template<>
    struct formatter<Physica::Real<Physica::Float64>, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        auto format(const Physica::Real<Physica::Float64>& obj, std::format_context& ctx) const {
            return std::format_to(ctx.out(), "{}", obj.toMachine());
        }
    };
}
