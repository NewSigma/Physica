/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Scalar/Rational.h"
#include <cassert>
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

using namespace Physica;

Rational::Rational(const Integer& i)
        : numerator(i), denominator(1) {}

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
    assert(!r.isZero());
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
    auto i = gcd<Integer, true>(numerator, denominator);
    numerator /= i;
    denominator /= i;
}

void Rational::swap(Rational& __restrict r) noexcept {
    assert(this != &r && "[Error]: Self swap is likely a bug");
    numerator.swap(r.numerator);
    denominator.swap(r.denominator);
}
