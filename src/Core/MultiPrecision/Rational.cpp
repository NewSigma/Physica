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
#include "Physica/Core/MultiPrecision/Rational.h"

namespace Physica::Core {
    Rational::Rational(const Integer& numerator_, const Integer& denominator_)
                         : numerator(numerator_), denominator(denominator_) {}

    Rational::Rational(const Rational& r)
                         : numerator(r.numerator), denominator(r.denominator) {}

    Rational::Rational(Rational&& r) noexcept
                         : numerator(std::move(r.numerator)), denominator(std::move(r.denominator)) {}

    Rational& Rational::operator=(const Rational& r) {
        if (this != &r) {
            numerator = r.numerator;
            denominator = r.denominator;
        }
        return *this;
    }

    Rational& Rational::operator=(Rational&& r) noexcept {
        numerator = std::move(r.numerator);
        denominator = std::move(r.denominator);
        return *this;
    }
}