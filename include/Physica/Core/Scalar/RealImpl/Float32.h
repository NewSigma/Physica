/*
 * Copyright 2024-2026 Weibo He.
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
#include <istream>
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
    class PHYSICA_API Real<Float32> : public ScalarBase<Real<Float32>>, public CRCoro<Real<Float32>> {
        using This = Real<Float32>;
        using Base = ScalarBase<This>;
    public:
        using device_obj_type = This;
    private:
        float f;
    public:
        constexpr Real() = default;
        [[gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr Real(std::floating_point auto f_) noexcept : f(f_) {}
        [[gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr Real(std::integral auto i) noexcept : f(i) {}
        Real(const Integer& i) : Real(float(double(i))) {}
        Real(const Rational& r) : Real(float(double(r))) {}
        template<Scalar T>
        __host__ __device__ explicit(Float32 < T::Prec) Real(const T& x) requires(!T::isComplex() && !Diffable<T>);
        constexpr Real(const Real&) = default;
        constexpr Real(Real&&) noexcept = default;
        constexpr ~Real() = default;
        /* Operators */
        using Base::operator=;
        using Base::operator<=>;
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
        [[nodiscard]] __host__ __device__ constexpr bool operator==(const Real& other) const noexcept { return f == other.f; }
        [[nodiscard]] __host__ __device__ constexpr auto operator<=>(const Real& other) const noexcept;
        [[nodiscard]] __host__ __device__ explicit operator float() const noexcept { return f; }
        [[nodiscard]] __host__ __device__ explicit operator double() const noexcept { return f; }
        [[nodiscard]] __host__ __device__ constexpr Real operator+(const Real& x) const noexcept;
        [[nodiscard]] __host__ __device__ constexpr Real operator-(const Real& x) const noexcept;
        [[nodiscard]] __host__ __device__ constexpr Real operator*(const Real& x) const noexcept;
        [[nodiscard]] __host__ __device__ constexpr Real operator/(const Real& x) const noexcept;
        [[nodiscard]] __host__ __device__ Real operator<<(int i) const { return Real(std::ldexp(f, i)); }
        [[nodiscard]] __host__ __device__ Real operator>>(int i) const { return Real(std::ldexp(f, -i)); }
        [[nodiscard]] __host__ __device__ Real operator-() const noexcept { return Real(-f); }
        friend std::istream& operator>>(std::istream& is, Real& scalar);
        /* Operations */
        [[nodiscard]] __host__ __device__ inline Real mod() const noexcept;
        [[nodiscard]] __host__ __device__ Real stripSignificand() const noexcept;
        void dump() const noexcept;

        using Base::random_uniform;
        using Base::random_normal;
        __host__ __device__ void swap(Real& __restrict x) noexcept { std::swap(f, x.f); }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr float toMachine() const noexcept { return f; }
        [[nodiscard]] __host__ __device__ constexpr float toMKL() const noexcept { return toMachine(); }
        [[nodiscard]] __host__ __device__ constexpr float toCUDA() const noexcept { return toMachine(); }
        [[nodiscard]] __host__ __device__ constexpr bool isZero() const noexcept { return f == 0; }
        [[nodiscard]] __host__ __device__ constexpr bool isSubNormal() const noexcept;
        [[nodiscard]] __host__ __device__ bool isPositive() const noexcept { return f > 0; }
        [[nodiscard]] __host__ __device__ bool isNegative() const noexcept { return f < 0; }
        [[nodiscard]] __host__ __device__ inline bool isFinite() const noexcept;
        [[nodiscard]] __host__ __device__ inline bool isInfinity() const noexcept;
        /* Static Members */
        [[nodiscard]] constexpr static Real nan() noexcept;
        template<RNG R>
        [[nodiscard]] static Real random_uniform() noexcept;
        template<RNG R>
        [[nodiscard]] static Real random_normal() noexcept;
        template<RNG R>
        [[nodiscard]] static Real random_normal(GaussRandomPool<This, R>& pool) noexcept { return pool(); }
        template<RNG R>
        [[nodiscard]] static Real random_any(auto& distribution) noexcept { return Real(distribution(R::getInstance())); }
    #ifdef PHYSICA_HDF5
        [[nodiscard]] static H5Type dtype_hdf5() noexcept { return H5Type::get<float>(); }
    #endif
    #ifdef PHYSICA_MPI
        [[nodiscard]] static MPI_Datatype dtype_mpi() noexcept { return MPI_FLOAT; }
    #endif
    };

    template<Scalar T>
    __host__ __device__ Real<Float32>::Real(const T& x) requires(!T::isComplex() && !Diffable<T>) : f(float(x)) {}

    __host__ __device__ constexpr auto Real<Float32>::operator<=>(const Real& other) const noexcept {
        auto order = f <=> other.f;
        if (order == std::partial_ordering::equivalent)
            return std::strong_ordering::equivalent;
        if (order == std::partial_ordering::less)
            return std::strong_ordering::less;
        if (order == std::partial_ordering::greater)
            return std::strong_ordering::greater;
        unreachable("Encounter NAN");
    }

    __host__ __device__ constexpr auto Real<Float32>::operator+(const Real& x) const noexcept -> Real {
        return Real(f + x.f);
    }

    __host__ __device__ constexpr auto Real<Float32>::operator-(const Real& x) const noexcept -> Real {
        return Real(f - x.f);
    }

    __host__ __device__ constexpr auto Real<Float32>::operator*(const Real& x) const noexcept -> Real {
        return Real(f * x.f);
    }

    __host__ __device__ constexpr auto Real<Float32>::operator/(const Real& x) const noexcept -> Real {
        assert((!x.isSubNormal() || x.isInfinity()) && "[Error]: Division overflow");
        return Real(f / x.f);
    }

    __host__ __device__ inline Real<Float32> Real<Float32>::mod() const noexcept {
        float buffer{};
        return std::modf(toMachine(), &buffer);
    }

    [[clang::no_sanitize("numerical")]] __host__ __device__ constexpr bool Real<Float32>::isSubNormal() const noexcept {
        return !__builtin_isnormal(f); // Use builtin to help no_sanitize
    }

    __host__ __device__ inline bool Real<Float32>::isFinite() const noexcept {
        return std::isfinite(f);
    }

    __host__ __device__ inline bool Real<Float32>::isInfinity() const noexcept {
        return std::isinf(f);
    }

    constexpr Real<Float32> Real<Float32>::nan() noexcept {
        return std::numeric_limits<float>::quiet_NaN();
    }

    template<RNG R>
    auto Real<Float32>::random_uniform() noexcept -> This {
        return Real(std::generate_canonical<float, std::numeric_limits<float>::digits>(R::getInstance()));
    }

    template<RNG R>
    auto Real<Float32>::random_normal() noexcept -> This {
        std::normal_distribution<float> dist{};
        return Real(dist(R::getInstance()));
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
    struct PHYSICA_API formatter<Physica::Real<Physica::Float32>, char> {
        constexpr static auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
        static auto format(Physica::Real<Physica::Float32> obj, auto& ctx) {
            if (obj.isZero())
                return std::format_to(ctx.out(), "0");
            return std::format_to(ctx.out(), "{}", obj.toMachine());
        }
    };
}

namespace Physica {
    inline std::ostream& operator<<(std::ostream& os, Real<Float32> x) {
        return os << std::format("{}", x.toMachine());
    }

    inline std::istream& operator>>(std::istream& is, Real<Float32>& scalar) {
        is >> scalar.f;
        return is;
    }
}
