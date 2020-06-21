/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SCALAR_H
#define PHYSICA_SCALAR_H

#include <iostream>
#include <cstring>
#include <cmath>
#include <QtCore/qlogging.h>
#include "Core/Header/SystemBits.h"
#include "Core/Header/Utils/MetaSupport.h"
#include "Core/Header/Const.h"
#include "Core/Header/MultiPrecition/ElementaryFunction.h"

namespace Physica::Core {
    /*!
     * \Scalar is a advanced float type that supports multiple precision and error track,
     * which is also compatible with float and double.
     */
    template<size_t maxPrecision = GlobalPrecision, bool errorTrack = true> class Scalar;
    //TODO operator+-*/ for different Scalars.
    template<size_t maxPrecision>
    class Scalar<maxPrecision, false> {
    protected:
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
        template<size_t maxPrecision2>
        Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> operator+(const Scalar<maxPrecision2, false>& s) const;
        template<size_t maxPrecision2>
        Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> operator-(const Scalar<maxPrecision2, false>& s) const;
        template<size_t maxPrecision2>
        Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> operator*(const Scalar<maxPrecision2, false>& s) const;
        template<size_t maxPrecision2>
        Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> operator/(const Scalar<maxPrecision2, false>& s) const;
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
        template<size_t maxPrecision2>
        inline static Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> add (
                const Scalar& s1, const Scalar<maxPrecision2, false>& s2);
        template<size_t maxPrecision2>
        inline static Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> sub (
                const Scalar& s1, const Scalar<maxPrecision2, false>& s2);
        template<size_t maxPrecision2>
        inline static Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> mul (
                const Scalar& s1, const Scalar<maxPrecision2, false>& s2);
        template<size_t maxPrecision2>
        inline static Scalar<META_MAX<maxPrecision, maxPrecision2>::value, false> div (
                const Scalar& s1, const Scalar<maxPrecision2, false>& s2);
        inline static bool cutLength(Scalar& s);
        inline static void cutZero(Scalar& s);
    };

    template<size_t maxPrecision>
    class Scalar<maxPrecision, true> : public Scalar<maxPrecision, false> {
        using Scalar<maxPrecision, false>::byte;
        using Scalar<maxPrecision, false>::length;
        using Scalar<maxPrecision, false>::power;
        //Accuracy
        ScalarUnit a;
    public:
        Scalar() noexcept;
        Scalar(int length, int power, ScalarUnit a = 0) noexcept;
        Scalar(const Scalar& s);
        Scalar(Scalar&& s) noexcept;
        explicit Scalar(SignedScalarUnit unit, ScalarUnit a = 0);
        explicit Scalar(double d, ScalarUnit a = 0);
        explicit Scalar(const char* s, ScalarUnit a = 0);
        explicit Scalar(const wchar_t* s, ScalarUnit a = 0);
        /* Operators */
        Scalar& operator=(const Scalar& s);
        Scalar& operator=(Scalar&& s) noexcept;
        template<size_t maxPrecision2>
        Scalar<META_MAX<maxPrecision, maxPrecision2>::value, true> operator+(const Scalar<maxPrecision2, true>& s) const;
        template<size_t maxPrecision2>
        Scalar<META_MAX<maxPrecision, maxPrecision2>::value, true> operator-(const Scalar<maxPrecision2, true>& s) const;
        template<size_t maxPrecision2>
        Scalar<META_MAX<maxPrecision, maxPrecision2>::value, true> operator*(const Scalar<maxPrecision2, true>& s) const;
        template<size_t maxPrecision2>
        Scalar<META_MAX<maxPrecision, maxPrecision2>::value, true> operator/(const Scalar<maxPrecision2, true>& s) const;
        Scalar operator<<(int bits) const;
        Scalar operator>>(int bits) const;
        Scalar operator-() const;
        /* Helpers */
        Scalar& applyError(const Scalar& error);
        Scalar& toUnitA() noexcept { a = 1; return *this; }
        Scalar& clearA() noexcept { a = 0; return *this; }
        inline void swap(Scalar& s) noexcept;
        /* Getters */
        [[nodiscard]] ScalarUnit getA() const noexcept { return a; }
        [[nodiscard]] Scalar<maxPrecision, false> getAccuracy() const;
        [[nodiscard]] inline Scalar getMaximum() const;
        [[nodiscard]] inline Scalar getMinimum() const;
    };
    /////////////////////////////////////////////Float////////////////////////////////////////////////
    template<>
    class Scalar<0, false> {
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
        /* Helpers */
        void swap(Scalar& s) noexcept { std::swap(f, s.f); }
    };

    template<>
    class Scalar<0, true> : public Scalar<0, false> {
        using Scalar<0, false>::f;
        float a;
    public:
        inline Scalar();
        inline explicit Scalar(float f, float a = 0);
        /* Operators */
        /* Helpers */
        void swap(Scalar& s) noexcept;
        Scalar operator+(const Scalar& s) const { return Scalar(f + s.f, a + s.a); }
        Scalar operator-(const Scalar& s) const { return Scalar(f - s.f, a + s.a); }
        Scalar operator*(const Scalar& s) const { return Scalar(f * s.f, f * s.a + s.f * a + a * s.a); }
        Scalar operator/(const Scalar& s) const { return Scalar(f / s.f, fabsf((f * a + s.f * s.a) / (s.f * (s.f - s.a)))); }
        /* Getters */
        [[nodiscard]] float getA() const noexcept { return a; }
    };
    /////////////////////////////////////////////Double////////////////////////////////////////////////
    template<>
    class Scalar<1, false> {
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
        /* Helpers */
        void swap(Scalar& s) noexcept { std::swap(d, s.d); }
    };

    template<>
    class Scalar<1, true> : public Scalar<1, false> {
        using Scalar<1, false>::d;
        double a;
    public:
        inline Scalar();
        inline explicit Scalar(double d, double a = 0);
        /* Operators */
        /* Helpers */
        void swap(Scalar& s) noexcept;
        Scalar operator+(const Scalar& s) const { return Scalar(d + s.d, a + s.a); }
        Scalar operator-(const Scalar& s) const { return Scalar(d - s.d, a + s.a); }
        Scalar operator*(const Scalar& s) const { return Scalar(d * s.d, d * s.a + s.d * a + a * s.a); }
        Scalar operator/(const Scalar& s) const { return Scalar(d / s.d, fabs((d * a + s.d * s.a) / (s.d * (s.d - s.a)))); }
        /* Getters */
        [[nodiscard]] double getA() const noexcept { return a; }
    };
    /* typedefs for convenience */
    [[maybe_unused]] typedef Scalar<0> FloatScalar;
    [[maybe_unused]] typedef Scalar<1> DoubleScalar;
    [[maybe_unused]] typedef Scalar<> MultiScalar;
}

#include "ScalarImpl/MultiScalar.h"
#include "ScalarImpl/BasicCalc.h"

#endif
