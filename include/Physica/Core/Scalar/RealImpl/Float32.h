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
#include <iomanip>
#include "Physica/CRCoro.h"
#include "../Real.h"

namespace Physica {
    template<>
    class Traits<float32> {
    public:
        constexpr static FloatPrec Prec = Float32;
        constexpr static int Order = 0;
        constexpr static bool isComplex = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;

        using ValueType = float32;
        using ScalarType = ValueType;
        using PtrTy = ScalarType*;
        using ConstPtrTy = const ScalarType*;
        using RefTy = ScalarType&;
        using ConstRefTy = const ScalarType&;
        using RealType = ScalarType;
        using ComplexType = Complex<ScalarType>;
        using MachineType = float;
    };
}

namespace Physica {
    template<>
    class Real<Float32> : public ScalarBase<Real<Float32>>, public CRCoro<Real<Float32>> {
        using This = Real<Float32>;
        using Base = ScalarBase<This>;
    public:
        using device_obj_type = This;
    private:
        float f;
    public:
        constexpr Real() = default;
        __host__ __device__ constexpr Real(float f_) noexcept : f(f_) {}
        Real(const Integer& i) : Real(float(double(i))) {}
        Real(const Rational& r) : Real(float(double(r))) {}
        template<Scalar T>
        __host__ __device__ explicit(Float32 < T::Prec) Real(const T& x) requires(!T::isComplex && !Diffable<T>);
        constexpr Real(const Real&) = default;
        constexpr Real(Real&&) noexcept = default;
        constexpr ~Real() = default;
        /* Operators */
        using Base::operator=;
        using Base::operator>;
        using Base::operator<;
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
        __host__ __device__ explicit operator float() const noexcept { return f; }
        __host__ __device__ explicit operator double() const noexcept { return f; }
        __host__ __device__ Real operator+(const Real& s) const noexcept { return Real(f + s.f); }
        __host__ __device__ Real operator-(const Real& s) const noexcept { return Real(f - s.f); }
        __host__ __device__ Real operator*(const Real& s) const noexcept { return Real(f * s.f); }
        __host__ __device__ Real operator/(const Real& s) const noexcept;
        __host__ __device__ Real operator<<(int i) const { return Real(std::ldexp(f, i)); }
        __host__ __device__ Real operator>>(int i) const { return Real(std::ldexp(f, -i)); }
        __host__ __device__ Real operator-() const noexcept { return Real(-f); }
        __host__ __device__ bool operator>(const Real& s) const noexcept { return f > s.f; }
        __host__ __device__ bool operator<(const Real& s) const noexcept { return f < s.f; }
        __host__ __device__ bool operator==(const Real& s) const noexcept { return f == s.f; }
        PHYSICA_API friend std::istream& operator>>(std::istream& is, Real& scalar);
        /* Operations */
        [[nodiscard]] __host__ __device__ inline Real mod() const noexcept;
        __host__ __device__ void swap(Real& __restrict s) noexcept { std::swap(f, s.f); }
        /* Getters */
        [[nodiscard]] __host__ __device__ float toMachine() const noexcept { return f; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return f == 0; }
        [[nodiscard]] __host__ __device__ bool isPositive() const noexcept { return f > 0; }
        [[nodiscard]] __host__ __device__ bool isNegative() const noexcept { return f < 0; }
        [[nodiscard]] __host__ __device__ inline bool isFinite() const noexcept;
        /* Static Members */
        [[nodiscard]] inline static Real nan() noexcept;
        template<RNG R>
        [[nodiscard]] static Real random_uniform();
        template<RNG R>
        [[nodiscard]] static Real random_normal();
        template<RNG R>
        [[nodiscard]] static Real random_normal(GaussRandomPool<This, R>& pool) { return pool(); }
        template<RNG R>
        [[nodiscard]] static Real random_any(auto& distribution) { return Real(distribution(R::getInstance())); }
        [[nodiscard]] static const H5::DataType& dtype_hdf5() noexcept { return H5::PredType::NATIVE_FLOAT; }
    #ifdef PHYSICA_MPI
        [[nodiscard]] static MPI_Datatype dtype_mpi() noexcept { return MPI_FLOAT; }
    #endif
    };

    template<Scalar T>
    __host__ __device__ Real<Float32>::Real(const T& x) requires(!T::isComplex && !Diffable<T>) : f(float(x)) {}

    __host__ __device__ inline Real<Float32> Real<Float32>::operator/(const This& s) const noexcept {
        assert(!s.isZero() && "[Error]: Divide by zero");
        return Real(f / s.f);
    }

    __host__ __device__ inline Real<Float32> Real<Float32>::mod() const noexcept {
        float buffer;
        return std::modf(toMachine(), &buffer);
    }

    __host__ __device__ inline bool Real<Float32>::isFinite() const noexcept {
        using namespace std;
        return isfinite(f);
    }

    inline Real<Float32> Real<Float32>::nan() noexcept {
        return std::nanf("");
    }

    template<RNG R>
    Real<Float32> Real<Float32>::random_uniform() {
        return Real(std::generate_canonical<float, std::numeric_limits<float>::digits>(R::getInstance()));
    }

    template<RNG R>
    Real<Float32> Real<Float32>::random_normal() {
        std::normal_distribution<float> dist{};
        return Real(dist(R::getInstance()));
    }

    inline std::ostream& operator<<(std::ostream& os, const Real<Float32>& x) {
        return os << std::format("{}", x.toMachine());
    }

    inline std::istream& operator>>(std::istream& is, Real<Float32>& scalar) {
        is >> scalar.f;
        return is;
    }
}

namespace std {
    template<>
    struct numeric_limits<Physica::Real<Physica::Float32>> : 
    #ifdef PHYSICA_CUDA
        public ::cuda::std::numeric_limits<float>
    #else
        public numeric_limits<float>
    #endif
    {};

    template<>
    struct formatter<Physica::Real<Physica::Float32>, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        auto format(const Physica::Real<Physica::Float32>& obj, std::format_context& ctx) const {
            return std::format_to(ctx.out(), "{}", obj.toMachine());
        }
    };
}
