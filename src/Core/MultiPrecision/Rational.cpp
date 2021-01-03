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
#include <cassert>
#include "Physica/Core/MultiPrecision/Rational.h"
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

namespace Physica::Core {
    Rational::Rational(const Integer& i) : numerator(i), denominator(1) {}
    
    Rational::Rational(const Integer& numerator_, const Integer& denominator_)
                         : numerator(numerator_), denominator(denominator_) {
        assert(denominator.isPositive());
    }

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

    Rational Rational::operator+(const Rational& r) const {
        Rational result(numerator * r.denominator + denominator * r.numerator, denominator * r.denominator);
        result.simplify();
        return result;
    }

    Rational Rational::operator-(const Rational& r) const {
        Rational result(numerator * r.denominator - denominator * r.numerator, denominator * r.denominator);
        result.simplify();
        return result;
    }
    /**
     * Optimize: (a/b) * (c/d) = (ac/bd), it is not clear whether simplify between a and d
     * (as well as between b and c) will improve the performance.
     */
    Rational Rational::operator*(const Rational& r) const {
        Rational result(numerator * r.numerator, denominator * r.denominator);
        result.simplify();
        return result;
    }
    /**
     * Optimize: refer to the Optimize above.
     */
    Rational Rational::operator/(const Rational& r) const {
        if (Q_UNLIKELY(r.isZero()))
            throw DivideByZeroException();
        Integer numerator_ = numerator * r.denominator;
        Integer denominator_ = denominator * r.numerator;
        numerator_.setSign(Integer::matchSign(numerator_, denominator_));
        denominator_.toAbs();
        Rational result(numerator_, denominator_);
        result.simplify();
        return result;
    }

    void Rational::simplify() {
        if (numerator.isZero()) {
            denominator = 1;
            return;
        }
        Integer gcd = GCD::run(numerator, denominator, GCD::Euclidean);
        numerator /= gcd;
        denominator /= gcd;
    }

    void Rational::swap(Rational& r) noexcept {
        numerator.swap(r.numerator);
        denominator.swap(r.denominator);
    }
}