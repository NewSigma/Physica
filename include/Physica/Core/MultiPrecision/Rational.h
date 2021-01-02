/*
 * Copyright 2021 WeiBo He.
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
#ifndef PHYSICA_RATIONAL_H
#define PHYSICA_RATIONAL_H

#include "Integer.h"

namespace Physica::Core {
    class Rational {
        Integer numerator;
        //Always possible integer.
        Integer denominator;
    public:
        Rational(const Integer& numerator_, const Integer& denominator_);
        Rational(const Rational& r);
        Rational(Rational&& r) noexcept;
        ~Rational() = default;
        /* Operators */
        Rational& operator=(const Rational& r);
        Rational& operator=(Rational&& r) noexcept;
        inline Rational operator+(const Rational& r) const;
        inline Rational operator-(const Rational& r) const;
        inline Rational operator*(const Rational& r) const;
        inline Rational operator/(const Rational& r) const;
        inline Rational operator-() const;
        /* Helpers */
        Rational& toOpposite() noexcept { numerator.toOpposite(); return *this; }
        Rational& toAbs() noexcept { numerator.toAbs(); return *this; }
        void swap(Rational& r) noexcept;
        /* Getters */
        [[nodiscard]] bool isZero() const { return numerator.isZero(); }
        [[nodiscard]] bool isPositive() const { return numerator.isPositive(); }
        [[nodiscard]] bool isNegative() const { return numerator.isNegative(); }
    };
}

#endif
