/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Integer.h"

namespace Physica {
    class PHYSICA_API Rational {
        Integer numerator;
        //Always positive integer.
        Integer denominator;
    public:
        Rational(const Integer& i);
        Rational(const Integer& numerator_, const Integer& denominator_);
        Rational(const Rational& r);
        Rational(Rational&& r) noexcept;
        ~Rational() = default;
        /* Operators */
        Rational& operator=(const Rational& r);
        Rational& operator=(Rational&& r) noexcept;
        Rational operator+(const Rational& r) const;
        Rational operator-(const Rational& r) const;
        Rational operator*(const Rational& r) const;
        Rational operator/(const Rational& r) const;
        inline Rational operator-() const;
        void operator+=(const Rational& r) { *this = *this + r; }
        void operator-=(const Rational& r) { *this = *this - r; }
        void operator*=(const Rational& r) { *this = *this * r; }
        void operator/=(const Rational& r) { *this = *this / r; }
        explicit operator double() const { return double(numerator) / double(denominator); }
        /* Operations */
        Rational& toOpposite() noexcept { numerator.toOpposite(); return *this; }
        Rational& toAbs() noexcept { numerator.toAbs(); return *this; }
        void simplify();
        void swap(Rational& __restrict r) noexcept;
        /* Getters */
        [[nodiscard]] const Integer& getNumerator() const noexcept { return numerator; }
        [[nodiscard]] const Integer& getDenominator() const noexcept { return denominator; }
        [[nodiscard]] bool isZero() const { return numerator.isZero(); }
        [[nodiscard]] bool isPositive() const { return numerator.isPositive(); }
        [[nodiscard]] bool isNegative() const { return numerator.isNegative(); }
    };

    inline Rational Rational::operator-() const {
        return Rational(-numerator, denominator);
    }

    inline void swap(Rational& r1, Rational& r2) noexcept {
        r1.swap(r2);
    }
}

namespace std {
    template<>
    inline void swap<Physica::Rational>(
            Physica::Rational& __restrict r1,
            Physica::Rational& __restrict r2) noexcept {
        r1.swap(r2);
    }
}
