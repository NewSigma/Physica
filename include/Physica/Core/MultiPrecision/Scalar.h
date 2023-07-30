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
#include "Rational.h"
#include "Physica/Utils/Template/CRTPBase.h"
#include "Physica/Utils/CUDA/device_obj.cuh"

namespace Physica::Core {
    //Forward declarations
    template<class AnyScalar> class ComplexScalar;

    template<ScalarOption option>
    __host__ __device__ Scalar<option> square(const Scalar<option>& s);

    template<ScalarOption option>
    __host__ __device__ Scalar<option> sqrt(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> ln(const Scalar<option>& s);

    namespace Internal {
        template<class T>
        class Traits;

        template<ScalarOption option_>
        class Traits<Scalar<option_>> {
            using Helper = typename std::conditional<option_ == Float, float, double>::type;
        public:
            using ScalarType = Scalar<option_>;
            using RealType = ScalarType;
            using ComplexType = ComplexScalar<ScalarType>;
            using TrivialType = typename std::conditional<option_ == MultiPrecision, Scalar<option_>, Helper>::type;
            static constexpr ScalarOption option = option_;
            static constexpr bool isComplex = false;
        };
        /**
         * This class return a type that can exactly represent the two input scalars.
         */
        template<class AnyScalar1, class AnyScalar2>
        class BinaryScalarOpReturnType {
            static constexpr ScalarOption option = Traits<AnyScalar1>::option > Traits<AnyScalar2>::option
                                                                              ? Traits<AnyScalar1>::option
                                                                              : Traits<AnyScalar2>::option;
            static constexpr bool isComplex = Traits<AnyScalar1>::isComplex || Traits<AnyScalar2>::isComplex;
        public:
            using Type = typename std::conditional<isComplex, ComplexScalar<Scalar<option>>, Scalar<option>>::type;
        };

        template<ScalarOption option>
        class AbstractScalar;

        template<>
        class AbstractScalar<MultiPrecision> {
            protected:
                //Store effective digits using little endian standard.
                MPUnit* __restrict byte;
                /*
                * Length of byte = abs(length).
                * sign of length and sign of Scalar are same. (when Scalar != 0)
                *
                * Warning: length can not equal to INT_MIN, or length will not return the correct answer.
                *
                * Optimize: use the end position of byte instead of length may improve performance.
                */
                int length;
                /*
                * Number = (x0 +- a * (2 ^ __WORDSIZE) ^ (1 - length)) * (2 ^ __WORDSIZE) ^ power
                *
                * FixIt: We have not considered overflow of power in our codes elsewhere.
                */
                int power;
            public:
                AbstractScalar() noexcept;
                AbstractScalar(int length_, int power_);
                AbstractScalar(const AbstractScalar& s);
                AbstractScalar(AbstractScalar<MultiPrecision>&& s) noexcept;
                AbstractScalar(int i);
                AbstractScalar(SignedMPUnit unit);
                AbstractScalar(double d);
                AbstractScalar(const Integer& i);
                AbstractScalar(const Rational& r);
                explicit AbstractScalar(const char* s);
                explicit AbstractScalar(const wchar_t* s);
                ~AbstractScalar();
                /* Operators */
                MPUnit operator[](unsigned int index) const { assert(index < static_cast<unsigned int>(getSize())); return byte[index]; }
                explicit operator double() const;
                AbstractScalar operator-() const;
                /* Helpers */
                void swap(AbstractScalar& s) noexcept;
                static inline bool matchSign(const AbstractScalar& s1, const AbstractScalar& s2);
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
            protected:
                /**
                 * Degigned for performance,
                 * this constructor should only be called by addNoError(), addWithError and etc.
                 *
                 * \param byte
                 * byte must be allocated by malloc()
                 */
                AbstractScalar(MPUnit* byte_, int length_, int power_) : byte(byte_), length(length_), power(power_) {}
                /* Operators */
                AbstractScalar& operator=(const AbstractScalar& s);
                AbstractScalar& operator=(AbstractScalar&& s) noexcept;
                /* Helpers */
                AbstractScalar& toOpposite() noexcept { length = -length; return *this; }
                AbstractScalar& toAbs() noexcept { length = getSize(); return *this; }
                inline void cutZero();
                /* Static members */
                inline static Scalar<MultiPrecision> addNoError(const AbstractScalar& s1, const AbstractScalar& s2);
                inline static Scalar<MultiPrecision> subNoError(const AbstractScalar& s1, const AbstractScalar& s2);
                inline static Scalar<MultiPrecision> mulNoError(const AbstractScalar& s1, const AbstractScalar& s2);
                inline static Scalar<MultiPrecision> divNoError(const AbstractScalar& s1, const AbstractScalar& s2);
                inline static bool cutLength(Scalar<MultiPrecision>& s);
                /* Friends */
                friend class Core::Integer;
                template<ScalarOption option>
                __host__ __device__ friend Scalar<option> Core::square(const Scalar<option>& s);
                template<ScalarOption option>
                __host__ __device__ friend Scalar<option> Core::sqrt(const Scalar<option>& s);
                template<ScalarOption option>
                friend Scalar<option> Core::ln(const Scalar<option>& s);
        };

        template<>
        class AbstractScalar<Float> {
        protected:
            float f;
        public:
            AbstractScalar() = default;
            __host__ __device__ AbstractScalar(float f_) : f(f_) {}
            AbstractScalar(const AbstractScalar& s) = default;
            AbstractScalar(const Integer& i) : AbstractScalar(float(double(i))) {}
            AbstractScalar(const Rational& r) : AbstractScalar(float(double(r))) {}
            ~AbstractScalar() = default;
            /* Operators */
            __host__ __device__ explicit operator float() const { return f; }
            __host__ __device__ explicit operator double() const { return f; }
            friend std::istream& operator>>(std::istream& is, AbstractScalar& scalar);
            /* Helpers */
            static inline bool matchSign(const AbstractScalar& s1, const AbstractScalar& s2) { return (s1.f > 0 && s2.f > 0) || (s1.f < 0 && s2.f < 0); }
            /* Getters */
            [[nodiscard]] constexpr static ScalarOption getOption() { return Float; }
            [[nodiscard]] __host__ __device__ float getTrivial() const noexcept { return f; }
            [[nodiscard]] __host__ __device__ bool isZero() const noexcept { return f == 0; }
            [[nodiscard]] __host__ __device__ bool isPositive() const { return f > 0; }
            [[nodiscard]] __host__ __device__ bool isNegative() const { return f < 0; }
            [[nodiscard]] bool isInteger() const;
        protected:
            /* Helpers */
            AbstractScalar& toOpposite() noexcept { f = -f; return *this; }
            __host__ __device__ AbstractScalar& toAbs() noexcept { f = fabsf(f); return *this; }
            void swap(AbstractScalar& s) noexcept { std::swap(f, s.f); }
        };

        template<>
        class AbstractScalar<Double> {
        protected:
            double d;
        public:
            AbstractScalar() = default;
            __host__ __device__ AbstractScalar(double d_) : d(d_) {}
            AbstractScalar(const AbstractScalar& s) = default;
            AbstractScalar(const Integer& i) : AbstractScalar(double(i)) {}
            AbstractScalar(const Rational& r) : AbstractScalar(double(r)) {}
            ~AbstractScalar() = default;
            /* Operators */
            __host__ __device__ explicit operator float() const { return d; }
            __host__ __device__ explicit operator double() const { return d; }
            friend std::istream& operator>>(std::istream& is, AbstractScalar& scalar);
            /* Helpers */
            static inline bool matchSign(const AbstractScalar& s1, const AbstractScalar& s2) { return (s1.d > 0 && s2.d > 0) || (s1.d < 0 && s2.d < 0); }
            /* Getters */
            [[nodiscard]] constexpr static ScalarOption getOption() { return Double; }
            [[nodiscard]] __host__ __device__ double getTrivial() const noexcept { return d; }
            [[nodiscard]] __host__ __device__ bool isZero() const noexcept{ return d == 0; }
            [[nodiscard]] __host__ __device__ bool isPositive() const { return d > 0; }
            [[nodiscard]] __host__ __device__ bool isNegative() const { return d < 0; }
            [[nodiscard]] bool isInteger() const;
        protected:
            /* Helpers */
            AbstractScalar& toOpposite() noexcept { d = -d; return *this; }
            __host__ __device__ AbstractScalar& toAbs() noexcept { d = fabs(d); return *this; }
            void swap(AbstractScalar& s) noexcept { std::swap(d, s.d); }
        };
    }

    template<class Derived>
    class ScalarBase : public Utils::CRTPBase<Derived> {
    public:
        using ScalarType = Derived;
        using TrivialType = typename Internal::Traits<Derived>::TrivialType;
        using RealType = typename Internal::Traits<Derived>::RealType;
        using ComplexType = typename Internal::Traits<Derived>::ComplexType;
        static constexpr ScalarOption option = Internal::Traits<Derived>::option;
        static constexpr bool isComplex = Internal::Traits<Derived>::isComplex;

        [[nodiscard]] const RealType& getReal() const {
            if constexpr (isComplex)
                return this->getDerived().getReal();
            else
                return this->getDerived();
        }

        [[nodiscard]] RealType getImag() const {
            if constexpr (isComplex)
                return this->getDerived().getImag();
            else
                return Derived(0);
        }

        [[nodiscard]] ScalarType conjugate() const {
            if constexpr (isComplex)
                return this->getDerived().conjugate();
            else
                return getReal();
        }

        [[nodiscard]] ScalarType unit() const {
            if constexpr (isComplex)
                return this->getDerived().unit();
            else
                return ScalarType(getReal().isNegative() ? -1 : 1);
        }

        [[nodiscard]] RealType norm() const {
            return sqrt(squaredNorm());
        }

        [[nodiscard]] RealType squaredNorm() const {
            if constexpr (isComplex)
                return this->getDerived().squaredNorm();
            else
                return square(this->getDerived());
        }

        constexpr static int size() { return 1; }

        Derived& load(const TrivialType* p) {
            this->getDerived() = *p;
            return this->getDerived();
        }

        void store(TrivialType* p) const {
            *p = this->getDerived().getTrivial();
        }

        Derived& load_partial(int n, const TrivialType* p) {
            if (n)
                load(p);
            return this->getDerived();
        }

        void store_partial(int n, TrivialType* p) const {
            if (n)
                store(p);
        }

        void insert(int index, TrivialType value) {
            (void)index;
            this->getDerived() = ScalarType(value);
        }
    };

    template<class T>
    struct is_scalar : public std::is_base_of<ScalarBase<T>, T> {};

    template<class ScalarType> ScalarType relativeError(const ScalarType& scalar1, const ScalarType& scalar2);
    template<class ScalarType>
    bool scalarNear(const ScalarBase<ScalarType>& scalar1, const ScalarBase<ScalarType>& scalar2, double precision);
    template<ScalarOption option>
    std::ostream& operator<<(std::ostream& os, const Scalar<option>& s);
    template<ScalarOption option>
    inline Scalar<option> operator+(const Scalar<option>& s);
    template<ScalarOption option, class T>
    __host__ __device__ inline std::enable_if_t<std::is_convertible<T, Scalar<option>>::value, void> operator+=(Scalar<option>& s1, const T& s2);
    template<ScalarOption option, class T>
    __host__ __device__ inline std::enable_if_t<std::is_convertible<T, Scalar<option>>::value, void> operator-=(Scalar<option>& s1, const T& s2);
    template<ScalarOption option, class T>
    __host__ __device__ inline std::enable_if_t<std::is_convertible<T, Scalar<option>>::value, void> operator*=(Scalar<option>& s1, const T& s2);
    template<ScalarOption option, class T>
    __host__ __device__ inline std::enable_if_t<std::is_convertible<T, Scalar<option>>::value, void> operator/=(Scalar<option>& s1, const T& s2);
    template<ScalarOption option>
    inline void operator^=(Scalar<option>& s1, const Scalar<option>& s2);
    template<ScalarOption option, bool errorTrack>
    inline void operator<<=(Scalar<option>& s, int bits);
    template<ScalarOption option, bool errorTrack>
    inline void operator>>=(Scalar<option>& s, int bits);
    template<ScalarOption option>
    __host__ __device__ inline bool operator>=(const Internal::AbstractScalar<option>& s1, const Internal::AbstractScalar<option>& s2);
    template<ScalarOption option>
    __host__ __device__ inline bool operator<=(const Internal::AbstractScalar<option>& s1, const Internal::AbstractScalar<option>& s2);
    template<ScalarOption option>
    __host__ __device__ inline bool operator!= (const Internal::AbstractScalar<option>& s1, const Internal::AbstractScalar<option>& s2);
    template<ScalarOption option>
    inline void swap(Scalar<option>& s1, Scalar<option>& s2) noexcept;
    /////////////////////////////////////////////MultiPrecision////////////////////////////////////////////////
    template<>
    class Scalar<MultiPrecision> : public ScalarBase<Scalar<MultiPrecision>>, public Internal::AbstractScalar<MultiPrecision> {
    public:
        using Base = Internal::AbstractScalar<MultiPrecision>;
        using ScalarType = Scalar<MultiPrecision>;
    public:
        Scalar() : Base() {}
        Scalar(int length_, int power_) : Base(length_, power_) {}
        Scalar(int i) : Base(i) {}
        Scalar(SignedMPUnit unit) : Base(unit) {}
        Scalar(double d) : Base(d) {}
        Scalar(const Integer& i) : Base(i) {}
        Scalar(const Rational& r) : Base(r) {}
        explicit Scalar(const char* s) : Base(s) {}
        explicit Scalar(const wchar_t* s) : Base(s) {}
        Scalar(const Scalar& s) = default;
        Scalar(Scalar&& s) noexcept = default;
        ~Scalar() = default;
        /* Operators */
        Scalar& operator=(const Scalar& s) = default;
        Scalar& operator=(Scalar&& s) noexcept = default;
        Scalar operator+(const Scalar& s) const;
        Scalar operator-(const Scalar& s) const;
        Scalar operator*(const Scalar& s) const;
        Scalar operator/(const Scalar& s) const;
        Scalar operator<<(int bits) const;
        Scalar operator>>(int bits) const;
        Scalar operator-() const;
        /* Helpers */
        Scalar& toOpposite() noexcept { return static_cast<Scalar&>(Base::toOpposite()); }
        Scalar& toAbs() noexcept { return static_cast<Scalar&>(Base::toAbs()); }
    protected:
        Scalar(MPUnit* byte_, int length_, int power_) : AbstractScalar(byte_, length_, power_) {}
        /* Friends */
        friend class Internal::AbstractScalar<MultiPrecision>;
    };
    static_assert(sizeof(Scalar<MultiPrecision>) == sizeof(Internal::AbstractScalar<MultiPrecision>), "Algorithms are based on this assumption.");

    inline Scalar<MultiPrecision>& operator++(Scalar<MultiPrecision>& s);
    inline Scalar<MultiPrecision>& operator--(Scalar<MultiPrecision>& s);
    inline Scalar<MultiPrecision> operator++(Scalar<MultiPrecision>& s, int);
    inline Scalar<MultiPrecision> operator--(Scalar<MultiPrecision>& s, int);
    /* Compare */
    bool absCompare(const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2);
    bool operator> (const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2);
    bool operator< (const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2);
    bool operator== (const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2);
    //IDEA: Comparisons between Scalar<MultiPrecision, true> may consider their accuracy.
    /////////////////////////////////////////////Float////////////////////////////////////////////////
    template<>
    class Scalar<Float> : public ScalarBase<Scalar<Float>>, public Internal::AbstractScalar<Float> {
    public:
        using Base = Internal::AbstractScalar<Float>;
        using ScalarType = Scalar<Float>;
        using device_obj_type = device_obj<ScalarType>;
    public:
        Scalar() = default;
        __host__ __device__ Scalar(float f_) : Base(f_) {}
        Scalar(const Integer& i) : Base(i) {}
        Scalar(const Rational& r) : Base(r) {}
        __host__ __device__ inline Scalar(const Scalar<Double>& s);
        Scalar(const Scalar& s) = default;
        //~Scalar() = default; /* Dynamic parallelism of CUDA 12.1 does not recognize that PlainStruct is trivial */
        /* Operators */
        __host__ __device__ Scalar operator+(const Scalar& s) const { return Scalar(f + s.f); }
        __host__ __device__ Scalar operator-(const Scalar& s) const { return Scalar(f - s.f); }
        __host__ __device__ Scalar operator*(const Scalar& s) const { return Scalar(f * s.f); }
        __host__ __device__ Scalar operator/(const Scalar& s) const { return Scalar(f / s.f); }
        Scalar operator<<(int i) const { return Scalar(f * std::pow(2, i)); }
        Scalar operator>>(int i) const { return Scalar(f / std::pow(2, i)); }
        __host__ __device__ Scalar operator-() const noexcept { return Scalar(-f); }
        /* Helpers */
        Scalar& toOpposite() noexcept { return static_cast<Scalar&>(Base::toOpposite()); }
        __host__ __device__ Scalar& toAbs() noexcept { return static_cast<Scalar&>(Base::toAbs()); }
        void swap(Scalar& s) noexcept { Base::swap(s); }
        /* Static Members */
        template<class RandomGenerator>
        [[nodiscard]] static Scalar random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static Scalar random_normal(RandomGenerator& gen);
    };

    inline Scalar<Float>& operator++(Scalar<Float>& s);
    inline Scalar<Float>& operator--(Scalar<Float>& s);
    inline Scalar<Float> operator++(Scalar<Float>& s, int);
    inline Scalar<Float> operator--(Scalar<Float>& s, int);
    /* Compare */
    inline bool absCompare(const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2);
    __host__ __device__ inline bool operator> (const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2);
    __host__ __device__ inline bool operator< (const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2);
    __host__ __device__ inline bool operator== (const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2);
    /////////////////////////////////////////////Double////////////////////////////////////////////////
    template<>
    class Scalar<Double> : public ScalarBase<Scalar<Double>>, public Internal::AbstractScalar<Double> {
    public:
        using Base = Internal::AbstractScalar<Double>;
        using ScalarType = Scalar<Double>;
        using device_obj_type = device_obj<ScalarType>;
    public:
        Scalar() = default;
        __host__ __device__ Scalar(double d_) : Base(d_) {}
        Scalar(const Integer& i) : Base(i) {}
        Scalar(const Rational& r) : Base(r) {}
        __host__ __device__ inline Scalar(const Scalar<Float>& s);
        Scalar(const Scalar& s) = default;
        //~Scalar() = default; /* Dynamic parallelism of CUDA 12.1 does not recognize that PlainStruct is trivial */
        /* Operators */
        __host__ __device__ explicit operator double() const { return d; }
        __host__ __device__ Scalar operator+(const Scalar& s) const { return Scalar(d + s.d); }
        __host__ __device__ Scalar operator-(const Scalar& s) const { return Scalar(d - s.d); }
        __host__ __device__ Scalar operator*(const Scalar& s) const { return Scalar(d * s.d); }
        __host__ __device__ Scalar operator/(const Scalar& s) const { return Scalar(d / s.d); }
        Scalar operator<<(int i) const { return Scalar(d * std::pow(2, i)); }
        Scalar operator>>(int i) const { return Scalar(d / std::pow(2, i)); }
        __host__ __device__ Scalar operator-() const noexcept { return Scalar(-d); }
        /* Helpers */
        Scalar& toOpposite() noexcept { return static_cast<Scalar&>(Base::toOpposite()); }
        __host__ __device__ Scalar& toAbs() noexcept { return static_cast<Scalar&>(Base::toAbs()); }
        void swap(Scalar& s) noexcept { std::swap(d, s.d); }
        /* Static Members */
        template<class RandomGenerator>
        [[nodiscard]] static Scalar random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static Scalar random_normal(RandomGenerator& gen);
    };

    inline Scalar<Double>& operator++(Scalar<Double>& s);
    inline Scalar<Double>& operator--(Scalar<Double>& s);
    inline Scalar<Double> operator++(Scalar<Double>& s, int);
    inline Scalar<Double> operator--(Scalar<Double>& s, int);
    /* Compare */
    inline bool absCompare(const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2);
    __host__ __device__ inline bool operator> (const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2);
    __host__ __device__ inline bool operator< (const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2);
    __host__ __device__ inline bool operator== (const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2);

    template<ScalarOption option>
    inline Scalar<option> operator^(const Scalar<option>& s1, const Scalar<option>& s2);
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
    struct numeric_limits<Physica::Core::Scalar<Physica::Core::Float>> : public numeric_limits<float> {};

    template<>
    struct numeric_limits<Physica::Core::Scalar<Physica::Core::Double>> : public numeric_limits<double> {};
}

#include "ScalarImpl/ScalarImpl.h"
