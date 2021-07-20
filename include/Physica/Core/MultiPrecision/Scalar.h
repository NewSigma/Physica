/*
 * Copyright 2019-2021 WeiBo He.
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
#ifndef PHYSICA_SCALAR_H
#define PHYSICA_SCALAR_H

#include <cmath>
#include <qglobal.h>
#include <ostream>
#include "MultiPrecisionType.h"
#include "Rational.h"

namespace Physica::Core {
    //Forward declarations
    namespace Internal {
        template<ScalarType type> class ScalarAddSubExpressionHelper;
        template<ScalarType type, bool errorTrack> class ScalarAddSubExpression;
    }

    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> square(const Scalar<MultiPrecision, errorTrack>& s);

    template<ScalarType type>
    Scalar<type, false> sqrt(const Scalar<type, false>& s);

    template<ScalarType type>
    Scalar<type, true> sqrt(const Scalar<type, true>& s);

    template<ScalarType type>
    Scalar<type, false> ln(const Scalar<type, false>& s);

    template<ScalarType type>
    Scalar<type, true> ln(const Scalar<type, true>& s);

    namespace Internal {
        /**
         * This class return a type that can exactly represent the two input scalars.
         */
        template<class T1, class T2>
        class ScalarOperationReturnType;

        template<ScalarType type1, bool errorTrack1, ScalarType type2, bool errorTrack2>
        class ScalarOperationReturnType<Scalar<type1, errorTrack1>, Scalar<type2, errorTrack2>> {
            static constexpr ScalarType type = type1 > type2 ? type2 : type1;
            static constexpr bool errorTrack = errorTrack1 || errorTrack2;
        public:
            using Type = Scalar<type, errorTrack>;
        };

        template<ScalarType type>
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
                MPUnit operator[](unsigned int index) const { Q_ASSERT(index < static_cast<unsigned int>(getSize())); return byte[index]; }
                explicit operator double() const;
                AbstractScalar operator-() const;
                /* Helpers */
                void toInteger();
                void swap(AbstractScalar& s) noexcept;
                static inline bool matchSign(const AbstractScalar& s1, const AbstractScalar& s2);
                /* Getters */
                [[nodiscard]] constexpr static ScalarType getType() { return MultiPrecision; }
                [[nodiscard]] int getLength() const noexcept { return length; }
                [[nodiscard]] int getPower() const noexcept { return power; }
                [[nodiscard]] int getSize() const noexcept { return abs(length); }
                [[nodiscard]] bool isZero() const { return byte[getSize() - 1] == 0; }
                [[nodiscard]] bool isPositive() const { return !isZero() && length > 0; }
                [[nodiscard]] bool isNegative() const { return !isZero() && length < 0; }
                [[nodiscard]] bool isInteger() const { return getSize() - 1 == power; }
                /* Setters */
                void setPower(int i) noexcept { power = i; }
                void setByte(unsigned int index, MPUnit value) { Q_ASSERT(index < static_cast<unsigned int>(getSize())); byte[index] = value; }
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
                inline static Scalar<MultiPrecision, true> addWithError(const AbstractScalar& s1, const AbstractScalar& s2);
                inline static Scalar<MultiPrecision, false> addNoError(const AbstractScalar& s1, const AbstractScalar& s2);
                inline static Scalar<MultiPrecision, true> subWithError(const AbstractScalar& s1, const AbstractScalar& s2);
                inline static Scalar<MultiPrecision, false> subNoError(const AbstractScalar& s1, const AbstractScalar& s2);
                inline static Scalar<MultiPrecision, true> mulWithError(const AbstractScalar& s1, const AbstractScalar& s2);
                inline static Scalar<MultiPrecision, false> mulNoError(const AbstractScalar& s1, const AbstractScalar& s2);
                inline static Scalar<MultiPrecision, true> divWithError(const AbstractScalar& s1, const AbstractScalar& s2);
                inline static Scalar<MultiPrecision, false> divNoError(const AbstractScalar& s1, const AbstractScalar& s2);
                template<bool errorTrack>
                inline static bool cutLength(Scalar<MultiPrecision, errorTrack>& s);
                /* Friends */
                friend class Core::Integer;
                template<bool errorTrack>
                friend Scalar<MultiPrecision, errorTrack> Core::square(const Scalar<MultiPrecision, errorTrack>& s);
                template<ScalarType type>
                friend Scalar<type, false> Core::sqrt(const Scalar<type, false>& s);
                template<ScalarType type>
                friend Scalar<type, true> Core::sqrt(const Scalar<type, true>& s);
                template<ScalarType type>
                friend Scalar<type, false> Core::ln(const Scalar<type, false>& s);
                template<ScalarType type>
                friend Scalar<type, true> Core::ln(const Scalar<type, true>& s);
        };

        template<>
        class AbstractScalar<Float> {
        public:
            using TrivialType = float;
        protected:
            float f;
        public:
            AbstractScalar() : f(0) {}
            AbstractScalar(float f_) : f(f_) {}
            AbstractScalar(const AbstractScalar& s) = default;
            AbstractScalar(const Integer& i) : AbstractScalar(float(double(i))) {}
            AbstractScalar(const Rational& r) : AbstractScalar(float(double(r))) {}
            ~AbstractScalar() = default;
            /* Operators */
            explicit operator float() const { return f; }
            explicit operator double() const { return f; }
            /* Getters */
            [[nodiscard]] constexpr static ScalarType getType() { return Float; }
            [[nodiscard]] float getTrivial() const noexcept { return f; }
            [[nodiscard]] bool isZero() const { return f == 0; }
            [[nodiscard]] bool isPositive() const { return f > 0; }
            [[nodiscard]] bool isNegative() const { return f < 0; }
            [[nodiscard]] bool isInteger() const;
        protected:
            /* Helpers */
            AbstractScalar& toOpposite() noexcept { f = -f; return *this; }
            AbstractScalar& toAbs() noexcept { f = fabsf(f); return *this; }
            void toInteger() { modff(f, &f); }
            void swap(AbstractScalar& s) noexcept { std::swap(f, s.f); }
        };

        template<>
        class AbstractScalar<Double> {
        public:
            using TrivialType = double;
        protected:
            double d;
        public:
            AbstractScalar() : d(0) {}
            AbstractScalar(double d_) : d(d_) {}
            AbstractScalar(const AbstractScalar& s) = default;
            AbstractScalar(const Integer& i) : AbstractScalar(double(i)) {}
            AbstractScalar(const Rational& r) : AbstractScalar(double(r)) {}
            ~AbstractScalar() = default;
            /* Operators */
            explicit operator float() const { return d; }
            explicit operator double() const { return d; }
            /* Getters */
            [[nodiscard]] constexpr static ScalarType getType() { return Double; }
            [[nodiscard]] double getTrivial() const noexcept { return d; }
            [[nodiscard]] bool isZero() const { return d == 0; }
            [[nodiscard]] bool isPositive() const { return d > 0; }
            [[nodiscard]] bool isNegative() const { return d < 0; }
            [[nodiscard]] bool isInteger() const;
        protected:
            /* Helpers */
            AbstractScalar& toOpposite() noexcept { d = -d; return *this; }
            AbstractScalar& toAbs() noexcept { d = fabs(d); return *this; }
            void toInteger() { modf(d, &d); }
            void swap(AbstractScalar& s) noexcept { std::swap(d, s.d); }
        };
    }

    template<>
    class Scalar<MultiPrecision, false> : public Internal::AbstractScalar<MultiPrecision> {
        using Base = Internal::AbstractScalar<MultiPrecision>;
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
        Scalar(const Scalar<MultiPrecision, true>& s);
        Scalar(Scalar<MultiPrecision, true>&& s);
        Scalar(const Scalar& s) = default;
        Scalar(Scalar&& s) noexcept = default;
        ~Scalar() = default;
        /* Operators */
        Scalar& operator=(const Scalar& s) = default;
        Scalar& operator=(Scalar&& s) noexcept = default;
        Scalar& operator=(const Scalar<MultiPrecision, true>& s);
        Scalar& operator=(Scalar<MultiPrecision, true>&& s) noexcept;
        Scalar operator*(const Scalar& s) const;
        Scalar operator/(const Scalar& s) const;
        Scalar<MultiPrecision> operator*(const Scalar<MultiPrecision>& s) const;
        Scalar<MultiPrecision> operator/(const Scalar<MultiPrecision>& s) const;
        Scalar operator<<(int bits) const;
        Scalar operator>>(int bits) const;
        Scalar operator-() const;
        /* Helpers */
        Scalar& toOpposite() noexcept { return static_cast<Scalar&>(Base::toOpposite()); }
        Scalar& toAbs() noexcept { return static_cast<Scalar&>(Base::toAbs()); }
        static inline Scalar getZero() { return Scalar(static_cast<SignedMPUnit>(0)); }
        static inline Scalar getOne() { return Scalar(static_cast<SignedMPUnit>(1)); }
        static inline Scalar getTwo() { return Scalar(static_cast<SignedMPUnit>(2)); }
        /* Getters */
        [[nodiscard]] constexpr static bool getErrorTrack() { return false; }
        /* Setters */
        Scalar& toUnitA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
        Scalar& clearA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
    protected:
        Scalar(MPUnit* byte_, int length_, int power_) : AbstractScalar(byte_, length_, power_) {}
        /* Friends */
        friend class Internal::AbstractScalar<MultiPrecision>;
        friend class Internal::ScalarAddSubExpressionHelper<MultiPrecision>;
        friend class Scalar<MultiPrecision, true>;
    };
    static_assert(sizeof(Scalar<MultiPrecision, false>) == sizeof(Internal::AbstractScalar<MultiPrecision>), "Algorithms are based on this assumption.");

    template<>
    class Scalar<MultiPrecision, true> final : public Internal::AbstractScalar<MultiPrecision> {
        using Base = Internal::AbstractScalar<MultiPrecision>;
    protected:
        //Accuracy
        MPUnit a;
    public:
        Scalar() noexcept;
        Scalar(int length_, int power_, MPUnit a_ = 0) noexcept;
        Scalar(const Scalar& s);
        Scalar(Scalar&& s) noexcept;
        explicit Scalar(int i, MPUnit a_ = 0);
        explicit Scalar(SignedMPUnit unit, MPUnit a_ = 0);
        explicit Scalar(double d, MPUnit a_ = 0);
        explicit Scalar(const Integer& i, MPUnit a_ = 0) : Base(i), a(a_) {}
        explicit Scalar(const Rational& r, MPUnit a_ = 0) : Base(r), a(a_) {}
        explicit Scalar(const char* s, MPUnit a_ = 0);
        explicit Scalar(const wchar_t* s, MPUnit a_ = 0);
        Scalar(const Scalar<MultiPrecision, false>& s, MPUnit a_ = 0);
        Scalar(Scalar<MultiPrecision, false>&& s, MPUnit a_ = 0);
        ~Scalar() = default;
        /* Operators */
        Scalar& operator=(const Scalar& s) = default;
        Scalar& operator=(Scalar&& s) noexcept = default;
        Scalar& operator=(const Scalar<MultiPrecision, false>& s);
        Scalar& operator=(Scalar<MultiPrecision, false>&& s) noexcept;
        Scalar<MultiPrecision, true> operator*(const Scalar<MultiPrecision, false>& s) const;
        Scalar<MultiPrecision, true> operator/(const Scalar<MultiPrecision, false>& s) const;
        Scalar<MultiPrecision, true> operator*(const Scalar<MultiPrecision, true>& s) const;
        Scalar<MultiPrecision, true> operator/(const Scalar<MultiPrecision, true>& s) const;
        Scalar operator<<(int bits) const;
        Scalar operator>>(int bits) const;
        Scalar operator-() const;
        /* Helpers */
        Scalar& toOpposite() noexcept { return static_cast<Scalar&>(Base::toOpposite()); }
        Scalar& toAbs() noexcept { return static_cast<Scalar&>(Base::toAbs()); }
        void toInteger();
        void swap(Scalar& s) noexcept;
        static inline Scalar getZero() { return Scalar(static_cast<SignedMPUnit>(0)); }
        static inline Scalar getOne() { return Scalar(static_cast<SignedMPUnit>(1)); }
        static inline Scalar getTwo() { return Scalar(static_cast<SignedMPUnit>(2)); }
        /* Getters */
        [[nodiscard]] constexpr static bool getErrorTrack() { return true; }
        [[nodiscard]] MPUnit getA() const noexcept { return a; }
        [[nodiscard]] Scalar<MultiPrecision, false> getAccuracy() const;
        [[nodiscard]] inline Scalar<MultiPrecision, false> getMaximum() const;
        [[nodiscard]] inline Scalar<MultiPrecision, false> getMinimum() const;
        /* Setters */
        Scalar& toUnitA() noexcept { a = 1; return *this; }
        Scalar& clearA() noexcept { a = 0; return *this; }
    private:
        Scalar(MPUnit* byte_, int length_, int power_, MPUnit a_ = 0)
                : AbstractScalar<MultiPrecision>(byte_, length_, power_), a(a_) {}
        Scalar& applyError(const Scalar<MultiPrecision, false>& error);
        /* Friends */
        friend class Internal::AbstractScalar<MultiPrecision>;
        friend class Internal::ScalarAddSubExpressionHelper<MultiPrecision>;
        friend class Scalar<MultiPrecision, false>;
        template<ScalarType type>
        friend Scalar<type, true> sqrt(const Scalar<type, true>& s);
        template<ScalarType type>
        friend Scalar<type, true> ln(const Scalar<type, true>& s);
    };
    /* Compare */
    bool absCompare(const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2);
    bool operator> (const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2);
    bool operator< (const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2);
    bool operator== (const Internal::AbstractScalar<MultiPrecision>& s1, const Internal::AbstractScalar<MultiPrecision>& s2);
    //IDEA: Comparisons between Scalar<MultiPrecision, true> may consider their accuracy.
    /////////////////////////////////////////////Float////////////////////////////////////////////////
    template<>
    class Scalar<Float, false> final : public Internal::AbstractScalar<Float> {
        using Base = Internal::AbstractScalar<Float>;
    public:
        Scalar() : Base() {}
        Scalar(float f_) : Base(f_) {}
        Scalar(const Integer& i) : Base(i) {}
        Scalar(const Rational& r) : Base(r) {}
        inline Scalar(const Scalar<Float, true>& s);
        inline explicit Scalar(const Scalar<Double, false>& s);
        Scalar(const Scalar& s) = default;
        ~Scalar() = default;
        /* Operators */
        Scalar operator*(const Scalar& s) const { return Scalar(f * s.f); }
        Scalar operator/(const Scalar& s) const { return Scalar(f / s.f); }
        inline Scalar<Float, true> operator*(const Scalar<Float, true>& s) const;
        inline Scalar<Float, true> operator/(const Scalar<Float, true>& s) const;
        Scalar operator<<(int i) const { return Scalar(f * pow(2, i)); }
        Scalar operator>>(int i) const { return Scalar(f / pow(2, i)); }
        Scalar operator-() const noexcept { return Scalar(-f); }
        /* Helpers */
        Scalar& toOpposite() noexcept { return static_cast<Scalar&>(Base::toOpposite()); }
        Scalar& toAbs() noexcept { return static_cast<Scalar&>(Base::toAbs()); }
        inline void toInteger();
        void swap(Scalar& s) noexcept { Base::swap(s); }
        static inline Scalar getZero() { return Scalar(0); }
        static inline Scalar getOne() { return Scalar(1); }
        static inline Scalar getTwo() { return Scalar(2); }
        /* Getters */
        [[nodiscard]] constexpr static bool getErrorTrack() { return false; }
        [[nodiscard]] constexpr static float getA() { return 0; }
        /* Setters */
        static void setA(float value) { assert(value == 0); /* Nothing, for the convenience of implement templates */ }
        Scalar& toUnitA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
        Scalar& clearA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
        /* Friends */
        friend class Scalar<Float, true>;
    };

    template<>
    class Scalar<Float, true> final : public Internal::AbstractScalar<Float> {
        using Base = Internal::AbstractScalar<Float>;
        float a;
    public:
        Scalar() : Base(), a(0) {}
        explicit Scalar(float f_, float a_ = 0) : Base(f_), a(fabsf(a_)) {}
        explicit Scalar(const Integer& i, float a_ = 0) : Base(i), a(a_) {}
        explicit Scalar(const Rational& r, float a_ = 0) : Base(r), a(a_) {}
        Scalar(const Scalar<Float, false>& s) : Base(s), a(0) {}
        Scalar(const Scalar& s) = default;
        ~Scalar() = default;
        /* Operators */
        Scalar operator*(const Scalar<Float, false>& s) const { return Scalar(f * s.f, s.f * getA()); }
        Scalar operator/(const Scalar<Float, false>& s) const { return Scalar(a / s.f); }
        Scalar operator*(const Scalar& s) const { return Scalar(f * s.f, f * s.a + s.f * a + a * s.a); }
        Scalar operator/(const Scalar& s) const { return Scalar(f / s.f, (f * a + s.f * s.a) / (s.f * (s.f - s.a))); }
        Scalar operator<<(int i) const { const float power = pow(2, i); return Scalar(f * power, a * power); }
        Scalar operator>>(int i) const { const float power = pow(2, i); return Scalar(f / power, a / power); }
        Scalar operator-() const noexcept { return Scalar(-f, a); }
        /* Helpers */
        Scalar& toOpposite() noexcept { return static_cast<Scalar&>(Base::toOpposite()); }
        Scalar& toAbs() noexcept { return static_cast<Scalar&>(Base::toAbs()); }
        inline void toInteger();
        void swap(Scalar& s) noexcept;
        static inline Scalar getZero() { return Scalar(0); }
        static inline Scalar getOne() { return Scalar(1); }
        static inline Scalar getTwo() { return Scalar(2); }
        /* Getters */
        [[nodiscard]] constexpr static bool getErrorTrack() { return false; }
        [[nodiscard]] float getA() const noexcept { return a; }
        /* Setters */
        void setA(float f_) { a = f_; }
        Scalar& toUnitA() noexcept { a = 1; return *this; }
        Scalar& clearA() noexcept { a = 0; return *this; }
        /* Friends */
        friend class Scalar<Float, false>;
    };
    /* Compare */
    inline bool absCompare(const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2);
    inline bool operator> (const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2);
    inline bool operator< (const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2);
    inline bool operator== (const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2);
    /////////////////////////////////////////////Double////////////////////////////////////////////////
    template<>
    class Scalar<Double, false> final : public Internal::AbstractScalar<Double> {
        using Base = Internal::AbstractScalar<Double>;
    public:
        Scalar() : Base() {}
        Scalar(double d_) : Base(d_) {}
        Scalar(const Integer& i) : Base(i) {}
        Scalar(const Rational& r) : Base(r) {}
        inline Scalar(const Scalar<Float, false>& s);
        Scalar(const Scalar<Double, true>& s);
        Scalar(const Scalar& s) = default;
        ~Scalar() = default;
        /* Operators */
        explicit operator double() const { return d; }
        Scalar operator*(const Scalar& s) const { return Scalar(d * s.d); }
        Scalar operator/(const Scalar& s) const { return Scalar(d / s.d); }
        inline Scalar<Double, true> operator*(const Scalar<Double, true>& s) const;
        inline Scalar<Double, true> operator/(const Scalar<Double, true>& s) const;
        Scalar operator<<(int i) const { return Scalar(d * pow(2, i)); }
        Scalar operator>>(int i) const { return Scalar(d / pow(2, i)); }
        Scalar operator-() const noexcept { return Scalar(-d); }
        /* Helpers */
        Scalar& toOpposite() noexcept { return static_cast<Scalar&>(Base::toOpposite()); }
        Scalar& toAbs() noexcept { return static_cast<Scalar&>(Base::toAbs()); }
        inline void toInteger();
        void swap(Scalar& s) noexcept { std::swap(d, s.d); }
        static inline Scalar getZero() { return Scalar(0); }
        static inline Scalar getOne() { return Scalar(1); }
        static inline Scalar getTwo() { return Scalar(2); }
        /* Getters */
        [[nodiscard]] constexpr static bool getErrorTrack() { return false; }
        [[nodiscard]] constexpr static double getA() { return 0; }
        /* Setters */
        static void setA(double value) { assert(value == 0); /* Nothing, for the convenience of implement templates */ }
        Scalar& toUnitA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
        Scalar& clearA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
        /* Friends */
        friend class Scalar<Double, true>;
    };

    template<>
    class Scalar<Double, true> final : public Internal::AbstractScalar<Double> {
        using Base = Internal::AbstractScalar<Double>;
        double a;
    public:
        Scalar() : Base(), a(0) {}
        explicit Scalar(double d_, double a_ = 0) : Base(d_), a(fabs(a_)) {}
        explicit Scalar(const Integer& i, double a_ = 0) : Base(i), a(a_) {}
        explicit Scalar(const Rational& r, double a_ = 0) : Base(r), a(a_) {}
        inline Scalar(const Scalar<Double, false>& s);
        Scalar(const Scalar& s) = default;
        ~Scalar() = default;
        /* Operators */
        Scalar operator*(const Scalar<Double, false>& s) const { return Scalar(d * s.d, s.d * getA()); }
        Scalar operator/(const Scalar<Double, false>& s) const { return Scalar(a / s.d); }
        Scalar operator*(const Scalar& s) const { return Scalar(d * s.d, d * s.a + s.d * a + a * s.a); }
        Scalar operator/(const Scalar& s) const { return Scalar(d / s.d, (d * a + s.d * s.a) / (s.d * (s.d - s.a))); }
        Scalar operator<<(int i) const { const double power = pow(2, i); return Scalar(d * power, a * power); }
        Scalar operator>>(int i) const { const double power = pow(2, i); return Scalar(d / power, a / power); }
        Scalar operator-() const noexcept { return Scalar(-d, a); }
        /* Helpers */
        Scalar& toOpposite() noexcept { return static_cast<Scalar&>(Base::toOpposite()); }
        Scalar& toAbs() noexcept { return static_cast<Scalar&>(Base::toAbs()); }
        inline void toInteger();
        void swap(Scalar& s) noexcept;
        static inline Scalar getZero() { return Scalar(0); }
        static inline Scalar getOne() { return Scalar(1); }
        static inline Scalar getTwo() { return Scalar(2); }
        /* Getters */
        [[nodiscard]] constexpr static bool getErrorTrack() { return true; }
        [[nodiscard]] double getA() const noexcept { return a; }
        /* Setters */
        void setA(double d_) { a = d_; }
        Scalar& toUnitA() noexcept { a = 1; return *this; }
        Scalar& clearA() noexcept { a = 0; return *this; }
        /* Friends */
        friend class Scalar<Double, false>;
    };
    /* Compare */
    inline bool absCompare(const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2);
    inline bool operator> (const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2);
    inline bool operator< (const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2);
    inline bool operator== (const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2);
    /* Output */
    template<ScalarType type, bool errorTrack>
    std::ostream& operator<<(std::ostream& os, const Scalar<type, errorTrack>& s);
    /* typedefs for convenience */
    [[maybe_unused]] typedef Scalar<Float> FloatScalar;
    [[maybe_unused]] typedef Scalar<Double> DoubleScalar;
    [[maybe_unused]] typedef Scalar<MultiPrecision> MultiScalar;
}

#include "Const.h"
#include "ScalarImpl/ScalarImpl.h"

#endif
