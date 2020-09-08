/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.

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
#include "ScalarType.h"

namespace Physica::Core {
    /*!
     * \Scalar is a advanced float type that supports multiple precision and error track,
     * which is also compatible with float and double.
     */
    template<ScalarType type = MultiPrecision, bool errorTrack = true> class Scalar;

    template<>
    class Scalar<MultiPrecision, false> {
        //Store effective digits.
        ScalarUnit* __restrict byte;
        /*
         * Length of byte = abs(length).
         * sign of length and sign of Scalar are same. (when Scalar != 0)
         *
         * Optimize: use the end position of byte instead length may improve performance.
         */
        int length;
        /*
         * Number = (x0 +- a * (2 ^ __WORDSIZE) ^ (1 - length)) * (2 ^ __WORDSIZE) ^ power
         * We have not considered overflow of power in our codes elsewhere.
         */
        int power;
    public:
        Scalar() noexcept;
        Scalar(int length, int power) noexcept;
        Scalar(const Scalar& s);
        Scalar(Scalar&& s) noexcept;
        Scalar(int i); //NOLINT Conversion is always available.
        Scalar(SignedScalarUnit unit); //NOLINT Conversion is always available.
        Scalar(double d); //NOLINT Conversion is always available.
        explicit Scalar(const char* s);
        explicit Scalar(const wchar_t* s);
        ~Scalar();
        /* Operators */
        Scalar& operator=(const Scalar& s);
        Scalar& operator=(Scalar&& s) noexcept;
        explicit operator double() const;
        ScalarUnit& operator[](unsigned int index) { return byte[index]; }
        const ScalarUnit& operator[](unsigned int index) const { return byte[index]; }
        Scalar operator+(const Scalar& s) const;
        Scalar operator-(const Scalar& s) const;
        Scalar operator*(const Scalar& s) const;
        Scalar operator/(const Scalar& s) const;
        Scalar<MultiPrecision> operator+(const Scalar<MultiPrecision>& s) const;
        Scalar<MultiPrecision> operator-(const Scalar<MultiPrecision>& s) const;
        Scalar<MultiPrecision> operator*(const Scalar<MultiPrecision>& s) const;
        Scalar<MultiPrecision> operator/(const Scalar<MultiPrecision>& s) const;
        Scalar operator<<(int bits) const;
        Scalar operator>>(int bits) const;
        Scalar operator-() const;
        /* Helpers */
        Scalar& toOpposite() noexcept { length = -length; return *this; }
        Scalar& toAbs() noexcept { length = getSize(); return *this; }
        void toInteger() noexcept;
        void swap(Scalar& s) noexcept;
        static inline bool matchSign(const Scalar& s1, const Scalar& s2);
        static inline Scalar getZero() { return Scalar(static_cast<SignedScalarUnit>(0)); }
        static inline Scalar getOne() { return Scalar(static_cast<SignedScalarUnit>(1)); }
        static inline Scalar getTwo() { return Scalar(static_cast<SignedScalarUnit>(2)); }
        /* Getters */
        [[nodiscard]] constexpr static ScalarType getType() { return MultiPrecision; }
        [[nodiscard]] constexpr static bool getErrorTrack() { return false; }
        [[nodiscard]] int getLength() const noexcept { return length; }
        [[nodiscard]] int getPower() const noexcept { return power; }
        [[nodiscard]] int getSize() const noexcept { return abs(length); }
        [[nodiscard]] bool isZero() const { return byte[getSize() - 1] == 0; }
        [[nodiscard]] bool isPositive() const { return !isZero() && length > 0; }
        [[nodiscard]] bool isNegative() const { return !isZero() && length < 0; }
        [[nodiscard]] bool isInteger() const { return getSize() == power + 1; }
        /* Setters */
        Scalar& toUnitA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
        Scalar& clearA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
    protected:
        //Should only be used in add(), sub(), mul() and div(). \byte must be allocated by malloc().
        Scalar(ScalarUnit* byte, int length, int power) : byte(byte), length(length), power(power) {}
        template<bool errorTrack>
        inline static Scalar<MultiPrecision, errorTrack> add (const Scalar& s1, const Scalar& s2);
        template<bool errorTrack>
        inline static Scalar<MultiPrecision, errorTrack> sub (const Scalar& s1, const Scalar& s2);
        template<bool errorTrack>
        inline static Scalar<MultiPrecision, errorTrack> mul (const Scalar& s1, const Scalar& s2);
        template<bool errorTrack>
        inline static Scalar<MultiPrecision, errorTrack> div (const Scalar& s1, const Scalar& s2);
        template<bool errorTrack>
        inline static bool cutLength(Scalar<MultiPrecision, errorTrack>& s);
        inline static void cutZero(Scalar& s);
        /* Friends */
        friend class Scalar<MultiPrecision, true>;
        template<bool errorTrack>
        friend Scalar<MultiPrecision, errorTrack> square(const Scalar<MultiPrecision, errorTrack>& s);
        template<ScalarType type>
        friend Scalar<type, false> sqrt(const Scalar<type, false>& s);
        template<ScalarType type>
        friend Scalar<type, true> sqrt(const Scalar<type, true>& s);
        template<ScalarType type>
        friend Scalar<type, false> ln(const Scalar<type, false>& s);
        template<ScalarType type>
        friend Scalar<type, true> ln(const Scalar<type, true>& s);
    };

    template<>
    class Scalar<MultiPrecision, true> : public Scalar<MultiPrecision, false> {
        //Accuracy
        ScalarUnit a;
    public:
        Scalar() noexcept;
        Scalar(int length, int power, ScalarUnit a = 0) noexcept;
        Scalar(const Scalar& s);
        Scalar(Scalar&& s) noexcept;
        Scalar(const Scalar<MultiPrecision, false>& s); //NOLINT Conversion is always available.
        Scalar(Scalar<MultiPrecision, false>&& s) noexcept; //NOLINT Conversion is always available.
        Scalar(int i, ScalarUnit a = 0); //NOLINT Conversion is always available.
        Scalar(SignedScalarUnit unit, ScalarUnit a = 0); //NOLINT Conversion is always available.
        Scalar(double d, ScalarUnit a = 0); //NOLINT Conversion is always available.
        explicit Scalar(const char* s, ScalarUnit a = 0);
        explicit Scalar(const wchar_t* s, ScalarUnit a = 0);
        /* Operators */
        Scalar& operator=(const Scalar& s);
        Scalar& operator=(Scalar&& s) noexcept;
        Scalar& operator=(const Scalar<MultiPrecision, false>& s);
        Scalar& operator=(Scalar<MultiPrecision, false>&& s) noexcept;
        Scalar<MultiPrecision, true> operator+(const Scalar<MultiPrecision, false>& s) const;
        Scalar<MultiPrecision, true> operator-(const Scalar<MultiPrecision, false>& s) const;
        Scalar<MultiPrecision, true> operator*(const Scalar<MultiPrecision, false>& s) const;
        Scalar<MultiPrecision, true> operator/(const Scalar<MultiPrecision, false>& s) const;
        Scalar<MultiPrecision, true> operator+(const Scalar<MultiPrecision, true>& s) const;
        Scalar<MultiPrecision, true> operator-(const Scalar<MultiPrecision, true>& s) const;
        Scalar<MultiPrecision, true> operator*(const Scalar<MultiPrecision, true>& s) const;
        Scalar<MultiPrecision, true> operator/(const Scalar<MultiPrecision, true>& s) const;
        Scalar operator<<(int bits) const;
        Scalar operator>>(int bits) const;
        Scalar operator-() const;
        /* Helpers */
        Scalar& applyError(const Scalar<MultiPrecision, false>& error);
        void toInteger() noexcept;
        void swap(Scalar& s) noexcept;
        static inline Scalar getZero() { return Scalar(static_cast<SignedScalarUnit>(0)); }
        static inline Scalar getOne() { return Scalar(static_cast<SignedScalarUnit>(1)); }
        static inline Scalar getTwo() { return Scalar(static_cast<SignedScalarUnit>(2)); }
        /* Getters */
        [[nodiscard]] constexpr static bool getErrorTrack() { return true; }
        [[nodiscard]] ScalarUnit getA() const noexcept { return a; }
        [[nodiscard]] Scalar<MultiPrecision, false> getAccuracy() const;
        [[nodiscard]] inline Scalar<MultiPrecision, false> getMaximum() const;
        [[nodiscard]] inline Scalar<MultiPrecision, false> getMinimum() const;
        /* Setters */
        Scalar& toUnitA() noexcept { a = 1; return *this; }
        Scalar& clearA() noexcept { a = 0; return *this; }
    private:
        //Should only be used in add(), sub(), mul() and div().
        Scalar(ScalarUnit* byte, int length, int power, ScalarUnit a = 0)
                : Scalar<MultiPrecision, false>(byte, length, power), a(a) {}
        /* Friends */
        friend class Scalar<MultiPrecision, false>;
    };
    /* Compare */
    bool absCompare(const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2);
    bool operator> (const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2);
    bool operator< (const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2);
    bool operator== (const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2);
    //IDEA: Comparisons between Scalar<MultiPrecision, true> may consider their accuracy.
    /////////////////////////////////////////////Float////////////////////////////////////////////////
    template<>
    class Scalar<Float, false> {
    protected:
        float f;
    public:
        inline Scalar();
        inline Scalar(float f); //NOLINT Intentional implicit conversions.
        Scalar(const Scalar& s) = default;
        /* Operators */
        explicit operator double() const { return f; }
        Scalar operator+(const Scalar& s) const { return Scalar(f + s.f); }
        Scalar operator-(const Scalar& s) const { return Scalar(f - s.f); }
        Scalar operator*(const Scalar& s) const { return Scalar(f * s.f); }
        Scalar operator/(const Scalar& s) const { return Scalar(f / s.f); }
        inline Scalar<Float, true> operator+(const Scalar<Float, true>& s) const;
        inline Scalar<Float, true> operator-(const Scalar<Float, true>& s) const;
        inline Scalar<Float, true> operator*(const Scalar<Float, true>& s) const;
        inline Scalar<Float, true> operator/(const Scalar<Float, true>& s) const;
        /* Helpers */
        Scalar& toOpposite() noexcept { f = -f; return *this; }
        Scalar& toAbs() noexcept { f = fabsf(f); return *this; }
        void swap(Scalar& s) noexcept { std::swap(f, s.f); }
        static inline Scalar getZero() { return Scalar(0); }
        static inline Scalar getOne() { return Scalar(1); }
        static inline Scalar getTwo() { return Scalar(2); }
        /* Getters */
        [[nodiscard]] constexpr static ScalarType getType() { return Float; }
        [[nodiscard]] constexpr static bool getErrorTrack() { return false; }
        [[nodiscard]] float getTrivial() const { return f; }
        [[nodiscard]] bool isZero() const { return f == 0; }
        [[nodiscard]] bool isPositive() const { return f > 0; }
        [[nodiscard]] bool isNegative() const { return f < 0; }
        [[nodiscard]] bool isInteger() const;
        /* Setters */
        static void setA(float b) { Q_UNUSED(b) /* Nothing, for the convenience of implement templates */ }
        Scalar& toUnitA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
        Scalar& clearA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
        /* Friends */
        friend class Scalar<Float, true>;
    };

    template<>
    class Scalar<Float, true> : public Scalar<Float, false> {
        float a;
    public:
        inline Scalar();
        inline explicit Scalar(float f, float a = 0);
        inline Scalar(const Scalar& s);
        /* Operators */
        Scalar operator+(const Scalar<Float, false>& s) const { return Scalar(f + s.f, getA()); }
        Scalar operator-(const Scalar<Float, false>& s) const { return Scalar(f - s.f, getA()); }
        Scalar operator*(const Scalar<Float, false>& s) const { return Scalar(f * s.f, s.f * getA()); }
        Scalar operator/(const Scalar<Float, false>& s) const { return Scalar(a / s.f); }
        Scalar operator+(const Scalar& s) const { return Scalar(f + s.f, a + s.a); }
        Scalar operator-(const Scalar& s) const { return Scalar(f - s.f, a + s.a); }
        Scalar operator*(const Scalar& s) const { return Scalar(f * s.f, f * s.a + s.f * a + a * s.a); }
        Scalar operator/(const Scalar& s) const { return Scalar(f / s.f, (f * a + s.f * s.a) / (s.f * (s.f - s.a))); }
        /* Helpers */
        void swap(Scalar& s) noexcept;
        static inline Scalar getZero() { return Scalar(0); }
        static inline Scalar getOne() { return Scalar(1); }
        static inline Scalar getTwo() { return Scalar(2); }
        /* Getters */
        [[nodiscard]] constexpr static bool getErrorTrack() { return false; }
        [[nodiscard]] float getA() const noexcept { return a; }
        /* Setters */
        void setA(float f) { a = f; }
        Scalar& toUnitA() noexcept { a = 1; return *this; }
        Scalar& clearA() noexcept { a = 0; return *this; }
        /* Friends */
        friend class Scalar<Float, false>;
    };
    /////////////////////////////////////////////Double////////////////////////////////////////////////
    template<>
    class Scalar<Double, false> {
    protected:
        double d;
    public:
        inline Scalar();
        inline Scalar(double d); //NOLINT Intentional implicit conversions.
        Scalar(const Scalar& s) = default;
        /* Operators */
        explicit operator double() const { return d; }
        Scalar operator+(const Scalar& s) const { return Scalar(d + s.d); }
        Scalar operator-(const Scalar& s) const { return Scalar(d - s.d); }
        Scalar operator*(const Scalar& s) const { return Scalar(d * s.d); }
        Scalar operator/(const Scalar& s) const { return Scalar(d / s.d); }
        inline Scalar<Double, true> operator+(const Scalar<Double, true>& s) const;
        inline Scalar<Double, true> operator-(const Scalar<Double, true>& s) const;
        inline Scalar<Double, true> operator*(const Scalar<Double, true>& s) const;
        inline Scalar<Double, true> operator/(const Scalar<Double, true>& s) const;
        /* Helpers */
        Scalar& toOpposite() noexcept { d = -d; return *this; }
        Scalar& toAbs() noexcept { d = fabs(d); return *this; }
        void swap(Scalar& s) noexcept { std::swap(d, s.d); }
        static inline Scalar getZero() { return Scalar(0); }
        static inline Scalar getOne() { return Scalar(1); }
        static inline Scalar getTwo() { return Scalar(2); }
        /* Getters */
        [[nodiscard]] constexpr static ScalarType getType() { return Double; }
        [[nodiscard]] constexpr static bool getErrorTrack() { return false; }
        [[nodiscard]] double getTrivial() const { return d; }
        [[nodiscard]] bool isZero() const { return d == 0; }
        [[nodiscard]] bool isPositive() const { return d > 0; }
        [[nodiscard]] bool isNegative() const { return d < 0; }
        [[nodiscard]] bool isInteger() const;
        /* Setters */
        static void setA(double d) { Q_UNUSED(d) /* Nothing, for the convenience of implement templates */ }
        Scalar& toUnitA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
        Scalar& clearA() noexcept { return *this; /* Nothing, for the convenience of implement templates */ }
        /* Friends */
        friend class Scalar<Double, true>;
    };

    template<>
    class Scalar<Double, true> : public Scalar<Double, false> {
        double a;
    public:
        inline Scalar();
        inline explicit Scalar(double d, double a = 0);
        inline Scalar(const Scalar& s);
        /* Operators */
        Scalar operator+(const Scalar<Double, false>& s) const { return Scalar(d + s.d, getA()); }
        Scalar operator-(const Scalar<Double, false>& s) const { return Scalar(d - s.d, getA()); }
        Scalar operator*(const Scalar<Double, false>& s) const { return Scalar(d * s.d, s.d * getA()); }
        Scalar operator/(const Scalar<Double, false>& s) const { return Scalar(a / s.d); }
        Scalar operator+(const Scalar& s) const { return Scalar(d + s.d, a + s.a); }
        Scalar operator-(const Scalar& s) const { return Scalar(d - s.d, a + s.a); }
        Scalar operator*(const Scalar& s) const { return Scalar(d * s.d, d * s.a + s.d * a + a * s.a); }
        Scalar operator/(const Scalar& s) const { return Scalar(d / s.d, (d * a + s.d * s.a) / (s.d * (s.d - s.a))); }
        /* Helpers */
        void swap(Scalar& s) noexcept;
        static inline Scalar getZero() { return Scalar(0); }
        static inline Scalar getOne() { return Scalar(1); }
        static inline Scalar getTwo() { return Scalar(2); }
        /* Getters */
        [[nodiscard]] constexpr static bool getErrorTrack() { return true; }
        [[nodiscard]] double getA() const noexcept { return a; }
        /* Setters */
        void setA(double d) { a = d; }
        Scalar& toUnitA() noexcept { a = 1; return *this; }
        Scalar& clearA() noexcept { a = 0; return *this; }
        /* Friends */
        friend class Scalar<Double, false>;
    };
    /* typedefs for convenience */
    [[maybe_unused]] typedef Scalar<Float> FloatScalar;
    [[maybe_unused]] typedef Scalar<Double> DoubleScalar;
    [[maybe_unused]] typedef Scalar<MultiPrecision> MultiScalar;
}

#include "Const.h"
#include "ScalarImpl/ScalarImpl.h"

#endif
