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
#pragma once

namespace Physica::Core {
    inline Rational Rational::operator+(const Rational& r) const {
        return Rational(numerator * r.denominator + denominator * r.numerator, denominator * r.denominator);
    }

    inline Rational Rational::operator-(const Rational& r) const {
        return Rational(numerator * r.denominator - denominator * r.numerator, denominator * r.denominator);
    }

    inline Rational Rational::operator*(const Rational& r) const {
        return Rational(numerator * r.numerator, denominator * r.denominator);
    }

    inline Rational Rational::operator/(const Rational& r) const {
        if (Q_UNLIKELY(r.isZero()))
            throw DivideByZeroException();
        return Rational(numerator * r.denominator, denominator * r.numerator);
    }

    inline Rational Rational::operator-() const {
        return Rational(-numerator, denominator);
    }

    inline void swap(Rational& r1, Rational& r2) noexcept {
        r1.swap(r2);
    }
}