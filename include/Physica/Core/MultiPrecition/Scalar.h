/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SCALAR_H
#define PHYSICA_SCALAR_H

#include <cmath>
#include "ScalarType.h"

namespace Physica::Core {
    enum ScalarType {
        Float,
        Double,
        MultiPrecision
    };
    /*!
     * \Scalar is a advanced float type that supports multiple precision and error track,
     * which is also compatible with float and double.
     */
    template<ScalarType type, bool errorTrack = true> class Scalar;

    template<>
    class Scalar<MultiPrecision, false> {
        //Store effective digits.
        ScalarUnit* __restrict byte;
        /*
         * Length of byte = abs(length).
         * sign of length and sign of Scalar are same. (when Scalar != 0)
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
        explicit Scalar(SignedScalarUnit unit);
        explicit Scalar(double d);
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
        void swap(Scalar& s) noexcept;
        /* Getters */
        [[nodiscard]] int getLength() const noexcept { return length; }
        [[nodiscard]] int getPower() const noexcept { return power; }
        [[nodiscard]] int getSize() const noexcept { return abs(length); }
        [[nodiscard]] bool isZero() const { return byte[getSize() - 1] == 0; }
        [[nodiscard]] bool isPositive() const { return !isZero() && length > 0; }
        [[nodiscard]] bool isNegative() const { return !isZero() && length < 0; }
        [[nodiscard]] bool isInteger() const { return getSize() == power + 1; }
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
        Scalar(const Scalar& s) = default;
        Scalar(Scalar&& s) noexcept;
        explicit Scalar(const Scalar<MultiPrecision, false>& s);
        explicit Scalar(Scalar<MultiPrecision, false>&& s) noexcept;
        explicit Scalar(SignedScalarUnit unit, ScalarUnit a = 0);
        explicit Scalar(double d, ScalarUnit a = 0);
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
        Scalar& toUnitA() noexcept { a = 1; return *this; }
        Scalar& clearA() noexcept { a = 0; return *this; }
        void swap(Scalar& s) noexcept;
        /* Getters */
        [[nodiscard]] ScalarUnit getA() const noexcept { return a; }
        [[nodiscard]] Scalar<MultiPrecision, false> getAccuracy() const;
        [[nodiscard]] inline Scalar<MultiPrecision, false> getMaximum() const;
        [[nodiscard]] inline Scalar<MultiPrecision, false> getMinimum() const;
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
        void swap(Scalar& s) noexcept { std::swap(f, s.f); }
        /* Getters */
        [[nodiscard]] float getTrivial() const { return f; }
        /* Friends */
        friend class Scalar<Float, true>;
    };

    template<>
    class Scalar<Float, true> : public Scalar<Float, false> {
        float a;
    public:
        inline Scalar();
        inline explicit Scalar(float f, float a = 0);
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
        /* Getters */
        [[nodiscard]] float getA() const noexcept { return a; }
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
        void swap(Scalar& s) noexcept { std::swap(d, s.d); }
        /* Getters */
        [[nodiscard]] double getTrivial() const { return d; }
        /* Friends */
        friend class Scalar<Double, true>;
    };

    template<>
    class Scalar<Double, true> : public Scalar<Double, false> {
        double a;
    public:
        inline Scalar();
        inline explicit Scalar(double d, double a = 0);
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
        /* Getters */
        [[nodiscard]] double getA() const noexcept { return a; }
        /* Friends */
        friend class Scalar<Double, false>;
    };
    /* typedefs for convenience */
    [[maybe_unused]] typedef Scalar<Float> FloatScalar;
    [[maybe_unused]] typedef Scalar<Double> DoubleScalar;
    [[maybe_unused]] typedef Scalar<MultiPrecision> MultiScalar;
}

#include "Physica/Core/Math/Const.h"
#include "Physica/Core/MultiPrecition/ScalarImpl/ScalarImpl.h"
#include "Physica/Core/MultiPrecition/ScalarImpl/BasicCalc.h"

#endif
