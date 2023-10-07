/*
 * Copyright 2019-2023 WeiBo He.
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

#include <cmath>
#include <ostream>
#include "MultiPrecisionType.h"
#ifdef PHYSICA_CUDA
    #include <cuda/std/limits>
#endif
#include "Rational.h"
#include "ScalarImpl/ScalarBase.h"
#include "Physica/Core/Exception/NotImplementedException.h"
#include "Physica/Core/IO/HDF5/HDF5.h"
#include "Physica/Utils/CUDA/device_obj.cuh"

namespace Physica::Core {
    //Forward declarations
    enum class DiffMode {
        Forward,
        Reverse
    };
    template<class ScalarType, DiffMode Mode> class Differentiable;
    template<class AnyScalar> class ComplexScalar;
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> abs(const Scalar<Option>& s);
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> square(const Scalar<Option>& s);
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> sqrt(const Scalar<Option>& s);
    template<ScalarOption Option> __host__ __device__ inline Scalar<Option> ln(const Scalar<Option>& s);

    namespace Internal {
        template<ScalarOption option_>
        class Traits<Scalar<option_>> {
            using Helper = typename std::conditional<option_ == Float, float, double>::type;
        public:
            using ScalarType = Scalar<option_>;
            using RealType = ScalarType;
            using ComplexType = ComplexScalar<ScalarType>;
            using TrivialType = typename std::conditional<option_ == MultiPrecision, Scalar<option_>, Helper>::type;
            using PlainScalar = ScalarType;
            static constexpr ScalarOption Option = option_;
            static constexpr bool isComplex = false;
            static constexpr bool isDifferentiable = false;
        };
        /**
         * This class return a type that can exactly represent the two input scalars.
         */
        template<class AnyScalar1, class AnyScalar2>
        class BinaryScalarOpReturnType {
            static_assert(Core::is_scalar<AnyScalar1>::value, "[Error]: This is not a scalar type");
            static_assert(Core::is_scalar<AnyScalar2>::value, "[Error]: This is not a scalar type");
            static_assert(!AnyScalar1::isComplex && !AnyScalar2::isComplex, "[Error]: This class applies to real scalar only");
            static_assert(!AnyScalar1::isDifferentiable && !AnyScalar2::isDifferentiable, "[Error]: This class applies to plain scalar only");
            static constexpr ScalarOption Option =
                    Traits<AnyScalar1>::Option > Traits<AnyScalar2>::Option ? Traits<AnyScalar1>::Option : Traits<AnyScalar2>::Option;
        public:
            using Type = Scalar<Option>;
        };
    }
    /////////////////////////////////////////////MultiPrecision////////////////////////////////////////////////
    template<>
    class Scalar<MultiPrecision> : public ScalarBase<Scalar<MultiPrecision>> {
    public:
        using ScalarType = Scalar<MultiPrecision>;
    private:
        //Store effective digits using little endian standard.
        MPUnit* __restrict byte;
        /**
         * Length of byte = abs(length).
         * sign of length and sign of Scalar are same. (when Scalar != 0)
         *
         * Warning: length can not equal to INT_MIN, or length will not return the correct answer.
         *
         * Optimize: use the end position of byte instead of length may improve performance.
         */
        int length;
        /**
         * Number = (x0 +- a * (2 ^ __WORDSIZE) ^ (1 - length)) * (2 ^ __WORDSIZE) ^ power
         *
         * FixIt: We have not considered overflow of power in our codes elsewhere.
         */
        int power;
    public:
        Scalar();
        Scalar(int length_, int power_);
        Scalar(int i);
        Scalar(SignedMPUnit unit);
        Scalar(double d);
        Scalar(const Integer& i);
        Scalar(const Rational& r);
        explicit Scalar(const char* s);
        explicit Scalar(const wchar_t* s);
        Scalar(const Scalar& s);
        Scalar(Scalar&& s) noexcept;
        ~Scalar();
        /* Operators */
        Scalar& operator=(Scalar s) noexcept;
        [[nodiscard]] MPUnit operator[](unsigned int index) const;
        [[nodiscard]] explicit operator double() const;
        [[nodiscard]] Scalar operator+(const Scalar& s) const;
        [[nodiscard]] Scalar operator-(const Scalar& s) const;
        [[nodiscard]] Scalar operator*(const Scalar& s) const;
        [[nodiscard]] Scalar operator/(const Scalar& s) const;
        [[nodiscard]] Scalar operator<<(int bits) const;
        [[nodiscard]] Scalar operator>>(int bits) const;
        [[nodiscard]] Scalar operator-() const;
        /* Operations */
        void swap(Scalar& obj) noexcept;
        /* Helpers */
        Scalar& toOpposite() noexcept { length = -length; return *this; }
        Scalar& toAbs() noexcept { length = getSize(); return *this; }
        /* Getters */
        [[nodiscard]] constexpr static ScalarOption getOption() { return MultiPrecision; }
        [[nodiscard]] int getLength() const noexcept { return length; }
        [[nodiscard]] int getPower() const noexcept { return power; }
        [[nodiscard]] int getSize() const noexcept { return std::abs(length); }
        [[nodiscard]] bool isZero() const { return byte[getSize() - 1] == 0; }
        [[nodiscard]] bool isPositive() const { return !isZero() && length > 0; }
        [[nodiscard]] bool isNegative() const { return !isZero() && length < 0; }
        [[nodiscard]] bool isInteger() const { return getSize() - 1 == power; }
        /* Setters */
        void setPower(int i) noexcept { power = i; }
        void setByte(unsigned int index, MPUnit value) { assert(index < static_cast<unsigned int>(getSize())); byte[index] = value; }
        /* Static members */
        static inline bool matchSign(const Scalar& s1, const Scalar& s2);
    private:
        /**
         * Degigned for performance,
         * this constructor should only be called by add(), sub() and etc.
         *
         * \param byte
         * byte must be allocated by malloc()
         */
        Scalar(MPUnit* byte_, int length_, int power_) : byte(byte_), length(length_), power(power_) {}
        /* Operations */
        inline void cutZero();
        /* Friends */
        friend class Integer;
        template<ScalarOption Option> __host__ __device__ friend Scalar<Option> square(const Scalar<Option>& s);
        template<ScalarOption Option> __host__ __device__ friend Scalar<Option> sqrt(const Scalar<Option>& s);
        template<ScalarOption Option> friend Scalar<Option> ln(const Scalar<Option>& s);
        /* Static members */
        inline static Scalar<MultiPrecision> add(const Scalar& s1, const Scalar& s2);
        inline static Scalar<MultiPrecision> sub(const Scalar& s1, const Scalar& s2);
        inline static Scalar<MultiPrecision> mul(const Scalar& s1, const Scalar& s2);
        inline static Scalar<MultiPrecision> div(const Scalar& s1, const Scalar& s2);
        inline static bool cutLength(Scalar<MultiPrecision>& s);
    };

    inline Scalar<MultiPrecision>& operator++(Scalar<MultiPrecision>& s);
    inline Scalar<MultiPrecision>& operator--(Scalar<MultiPrecision>& s);
    inline Scalar<MultiPrecision> operator++(Scalar<MultiPrecision>& s, int);
    inline Scalar<MultiPrecision> operator--(Scalar<MultiPrecision>& s, int);
    /* Compare */
    bool absCompare(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2);
    bool operator> (const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2);
    bool operator< (const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2);
    bool operator== (const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2);
    /////////////////////////////////////////////Float////////////////////////////////////////////////
    template<>
    class Scalar<Float> : public ScalarBase<Scalar<Float>> {
        using Base = ScalarBase<Scalar<Float>>;
    public:
        using ScalarType = Scalar<Float>;
        using device_obj_type = device_obj<ScalarType>;
    private:
        float f;
    public:
        Scalar() = default;
        __host__ __device__ Scalar(float f_) : f(f_) {}
        Scalar(const Integer& i) : Scalar(float(double(i))) {}
        Scalar(const Rational& r) : Scalar(float(double(r))) {}
        __host__ __device__ inline Scalar(const Scalar<Double>& s);
        Scalar(const Scalar&) = default;
        Scalar(Scalar&&) noexcept = default;
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
        Scalar operator<<(int i) const { return Scalar(f * std::pow(2, i)); }
        Scalar operator>>(int i) const { return Scalar(f / std::pow(2, i)); }
        __host__ __device__ Scalar operator-() const noexcept { return Scalar(-f); }
        __host__ __device__ bool operator>(const Scalar& s) const { return f > s.f; }
        __host__ __device__ bool operator<(const Scalar& s) const { return f < s.f; }
        __host__ __device__ bool operator==(const Scalar& s) const { return f == s.f; }
        friend std::istream& operator>>(std::istream& is, Scalar& scalar);
        /* Helpers */
        Scalar& toOpposite() noexcept { f = -f; return *this; }
        __host__ __device__ Scalar& toAbs() noexcept { *this = abs(*this); return *this; }
        void swap(Scalar& s) noexcept { std::swap(f, s.f); }
        /* Getters */
        [[nodiscard]] constexpr static ScalarOption getOption() { return Float; }
        [[nodiscard]] __host__ __device__ float getTrivial() const noexcept { return f; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return f == 0; }
        [[nodiscard]] __host__ __device__ bool isPositive() const { return f > 0; }
        [[nodiscard]] __host__ __device__ bool isNegative() const { return f < 0; }
        [[nodiscard]] bool isInteger() const;
        /* Static Members */
        static inline bool matchSign(const Scalar& s1, const Scalar& s2) { return (s1.f > 0 && s2.f > 0) || (s1.f < 0 && s2.f < 0); }
        template<class RandomGenerator>
        [[nodiscard]] static Scalar random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static Scalar random_normal(RandomGenerator& gen);
        [[nodiscard]] static const H5::DataType& getH5DataType() { return H5::PredType::NATIVE_FLOAT; }
    };

    inline Scalar<Float>& operator++(Scalar<Float>& s);
    inline Scalar<Float>& operator--(Scalar<Float>& s);
    inline Scalar<Float> operator++(Scalar<Float>& s, int);
    inline Scalar<Float> operator--(Scalar<Float>& s, int);
    /////////////////////////////////////////////Double////////////////////////////////////////////////
    template<>
    class Scalar<Double> : public ScalarBase<Scalar<Double>> {
        using Base = ScalarBase<Scalar<Double>>;
    public:
        using ScalarType = Scalar<Double>;
        using device_obj_type = device_obj<ScalarType>;
    private:
        double d;
    public:
        Scalar() = default;
        __host__ __device__ Scalar(double d_) : d(d_) {}
        Scalar(const Integer& i) : Scalar(double(i)) {}
        Scalar(const Rational& r) : Scalar(double(r)) {}
        __host__ __device__ inline Scalar(const Scalar<Float>& s);
        Scalar(const Scalar&) = default;
        Scalar(Scalar&&) noexcept = default;
        //~Scalar() = default; /* Dynamic parallelism of CUDA 12.1 does not recognize that PlainStruct is trivial */
        /* Operators */
        using Base::operator>;
        using Base::operator<;
        Scalar& operator=(const Scalar& obj) = default;
        Scalar& operator=(Scalar&& obj) noexcept = default;
        __host__ __device__ explicit operator float() const { return d; }
        __host__ __device__ explicit operator double() const { return d; }
        __host__ __device__ Scalar operator+(const Scalar& s) const { return Scalar(d + s.d); }
        __host__ __device__ Scalar operator-(const Scalar& s) const { return Scalar(d - s.d); }
        __host__ __device__ Scalar operator*(const Scalar& s) const { return Scalar(d * s.d); }
        __host__ __device__ Scalar operator/(const Scalar& s) const { return Scalar(d / s.d); }
        Scalar operator<<(int i) const { return Scalar(d * std::pow(2, i)); }
        Scalar operator>>(int i) const { return Scalar(d / std::pow(2, i)); }
        __host__ __device__ Scalar operator-() const noexcept { return Scalar(-d); }
        __host__ __device__ bool operator>(const Scalar& s) const { return d > s.d; }
        __host__ __device__ bool operator<(const Scalar& s) const { return d < s.d; }
        __host__ __device__ bool operator==(const Scalar& s) const { return d == s.d; }
        friend std::istream& operator>>(std::istream& is, Scalar& scalar);
        /* Helpers */
        Scalar& toOpposite() noexcept { d = -d; return *this; }
        __host__ __device__ Scalar& toAbs() noexcept { *this = abs(*this); return *this; }
        void swap(Scalar& s) noexcept { std::swap(d, s.d); }
        /* Getters */
        [[nodiscard]] constexpr static ScalarOption getOption() { return Double; }
        [[nodiscard]] __host__ __device__ double getTrivial() const noexcept { return d; }
        [[nodiscard]] __host__ __device__ bool isZero() const noexcept{ return d == 0; }
        [[nodiscard]] __host__ __device__ bool isPositive() const { return d > 0; }
        [[nodiscard]] __host__ __device__ bool isNegative() const { return d < 0; }
        [[nodiscard]] bool isInteger() const;
        /* Static Members */
        static inline bool matchSign(const Scalar& s1, const Scalar& s2) { return (s1.d > 0 && s2.d > 0) || (s1.d < 0 && s2.d < 0); }
        template<class RandomGenerator>
        [[nodiscard]] static Scalar random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static Scalar random_normal(RandomGenerator& gen);
        [[nodiscard]] static const H5::DataType& getH5DataType() { return H5::PredType::NATIVE_DOUBLE; }
    };

    inline Scalar<Double>& operator++(Scalar<Double>& s);
    inline Scalar<Double>& operator--(Scalar<Double>& s);
    inline Scalar<Double> operator++(Scalar<Double>& s, int);
    inline Scalar<Double> operator--(Scalar<Double>& s, int);

    template<ScalarOption Option>
    inline Scalar<Option> operator^(const Scalar<Option>& s1, const Scalar<Option>& s2);
    /* typedefs for convenience */
    [[maybe_unused]] typedef Scalar<Float> FloatScalar;
    [[maybe_unused]] typedef Scalar<Double> DoubleScalar;
    [[maybe_unused]] typedef Scalar<MultiPrecision> MultiScalar;
}

#include "Const.h"

namespace std {
    template<>
    struct numeric_limits<Physica::Core::Scalar<Physica::Core::MultiPrecision>> {
        using T = Physica::Core::Scalar<Physica::Core::MultiPrecision>;
    public:
        static T epsilon() noexcept {
            auto result = T(static_cast<SignedMPUnit>(1));
            result.setPower(1 - Physica::Core::GlobalPrecision);
            return result;
        }
    };

    template<>
    struct numeric_limits<Physica::Core::Scalar<Physica::Core::Float>> : 
    #ifdef PHYSICA_CUDA
        public ::cuda::std::numeric_limits<float>
    #else
        public numeric_limits<float>
    #endif
    {};

    template<>
    struct numeric_limits<Physica::Core::Scalar<Physica::Core::Double>> : 
    #ifdef PHYSICA_CUDA
        public ::cuda::std::numeric_limits<double>
    #else
        public numeric_limits<double>
    #endif
    {};
}

#include "ScalarImpl/ScalarImpl.h"
#include "SIMD.h"
