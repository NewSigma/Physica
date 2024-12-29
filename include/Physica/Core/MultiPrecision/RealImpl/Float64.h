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

#include <format>
#include "Physica/CRCoro.h"
#include "../Real.h"

namespace Physica {
    template<>
    class Traits<float64> {
    public:
        constexpr static ScalarOption Option = Float64;
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

namespace Physica::Core {
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
        __host__ __device__ explicit inline Real(const T& x);
        constexpr Real(const Real&) = default;
        constexpr Real(Real&&) noexcept = default;
        constexpr ~Real() = default;
        /* Operators */
        using Base::operator>;
        using Base::operator<;
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
        template<ScalarOption Op>
        __host__ __device__ This& operator=(Real<Op> x) noexcept;
        __host__ __device__ explicit operator float() const noexcept { return d; }
        __host__ __device__ explicit operator double() const noexcept { return d; }
        __host__ __device__ Real operator+(const Real& s) const noexcept { return Real(d + s.d); }
        __host__ __device__ Real operator-(const Real& s) const noexcept { return Real(d - s.d); }
        __host__ __device__ Real operator*(const Real& s) const noexcept { return Real(d * s.d); }
        __host__ __device__ Real operator/(const Real& s) const noexcept { return Real(d / s.d); }
        Real operator<<(int i) const { return Real(d * std::pow(2, i)); }
        Real operator>>(int i) const { return Real(d / std::pow(2, i)); }
        __host__ __device__ Real operator-() const noexcept { return Real(-d); }
        __host__ __device__ bool operator>(const Real& s) const noexcept { return d > s.d; }
        __host__ __device__ bool operator<(const Real& s) const noexcept { return d < s.d; }
        __host__ __device__ bool operator==(const Real& s) const noexcept { return d == s.d; }
        PHYSICA_API friend std::istream& operator>>(std::istream& is, Real& scalar);
        /* Operations */
        [[nodiscard]] inline Real mod() const noexcept;
        void swap(Real& __restrict s) noexcept { std::swap(d, s.d); }
        /* Getters */
        [[nodiscard]] constexpr static ScalarOption getOption() { return Float64; }
        [[nodiscard]] __host__ __device__ double toMachine() const noexcept { return d; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept{ return d == 0; }
        [[nodiscard]] __host__ __device__ bool isPositive() const noexcept { return d > 0; }
        [[nodiscard]] __host__ __device__ bool isNegative() const noexcept { return d < 0; }
        [[nodiscard]] __host__ __device__ inline bool isFinite() const noexcept;
        [[nodiscard]] bool isInteger() const;
        /* Static Members */
        [[nodiscard]] inline static Real nan() noexcept;
        template<RandomGenerator R>
        [[nodiscard]] inline static Real random_uniform();
        template<RandomGenerator R>
        [[nodiscard]] inline static Real random_normal();
        template<RandomGenerator R>
        [[nodiscard]] static Real random_normal(GaussRandomPool<This, R>& pool) { return pool(); }
        template<class Distribution, RandomGenerator R>
        [[nodiscard]] inline static Real random_any(Distribution& dist);
    #ifdef PHYSICA_HDF5
        [[nodiscard]] static const H5::DataType& getH5DataType() { return H5::PredType::NATIVE_DOUBLE; }
    #endif
    };

    template<Scalar T>
    __host__ __device__ inline Real<Float64>::Real(const T& x) : d(double(x)) {}

    template<ScalarOption Op>
    __host__ __device__ Real<Float64>& Real<Float64>::operator=(Real<Op> x) noexcept {
        return operator=(This(x));
    }

    inline Real<Float64> Real<Float64>::mod() const noexcept {
        double buffer;
        return std::modf(toMachine(), &buffer);
    }

    __host__ __device__ inline bool Real<Float64>::isFinite() const noexcept {
    #ifdef __CUDA_ARCH__
        return isfinite(d);
    #else
        return std::isfinite(d);
    #endif
    }

    inline Real<Float64> Real<Float64>::nan() noexcept {
        return std::nan("");
    }

    template<RandomGenerator R>
    inline Real<Float64> Real<Float64>::random_uniform() {
        std::uniform_real_distribution<double> dist{};
        return Real(dist(R::getInstance()));
    }

    template<RandomGenerator R>
    inline Real<Float64> Real<Float64>::random_normal() {
        std::normal_distribution<double> dist{};
        return Real(dist(R::getInstance()));
    }

    template<class Distribution, RandomGenerator R>
    inline Real<Float64> Real<Float64>::random_any(Distribution& dist) {
        return Real(dist(R::getInstance()));
    }

    inline std::ostream& operator<<(std::ostream& os, const Real<Float64>& s) {
        const auto lastPrec = os.precision();
        return os << std::setprecision(16) << double(s) << std::setprecision(lastPrec);
    }

    inline Real<Float64>& operator++(Real<Float64>& s) {
        s += float64(1.0);
        return s;
    }

    inline Real<Float64>& operator--(Real<Float64>& s) {
        s -= float64(1.0);
        return s;
    }

    inline Real<Float64> operator++(Real<Float64>& s, int) {
        float64 temp(s);
        s += float64(1.0);
        return temp;
    }

    inline Real<Float64> operator--(Real<Float64>& s, int) {
        float64 temp(s);
        s -= float64(1.0);
        return temp;
    }
}

namespace std {
    template<>
    struct numeric_limits<Physica::Core::Real<Physica::Core::Float64>> : 
    #ifdef PHYSICA_CUDA
        public ::cuda::std::numeric_limits<double>
    #else
        public numeric_limits<double>
    #endif
    {};

    template<>
    struct formatter<Physica::Core::Real<Physica::Core::Float64>, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        auto format(const Physica::Core::Real<Physica::Core::Float64>& obj, std::format_context& ctx) const {
            return std::format_to(ctx.out(), "{:.16G}", obj.toMachine());
        }
    };
}
