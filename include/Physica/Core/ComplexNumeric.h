/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_COMPLEXNUMERIC_H
#define PHYSICA_COMPLEXNUMERIC_H

#include "Scalar.h"

namespace Physica::Core {
    class ComplexNumeric : private Scalar {
        Scalar imagine;
    public:
        ComplexNumeric() = default;
        ComplexNumeric(Scalar n1, Scalar n2) noexcept;
        ComplexNumeric(const ComplexNumeric& c) = default;
        ComplexNumeric(ComplexNumeric&& c) noexcept;
        /* Operators */
        friend std::ostream& operator<<(std::ostream& os, const ComplexNumeric& n);
        ComplexNumeric& operator=(const ComplexNumeric& c);
        ComplexNumeric& operator=(ComplexNumeric&& c) noexcept;
        /* Helpers */
        void swap(ComplexNumeric& c) noexcept;
        /* Getters */
        [[nodiscard]] const Scalar& getReal() const { return *this; }
        [[nodiscard]] const Scalar& getImagine() const { return imagine; }
    };
    [[nodiscard]] Scalar norm(const ComplexNumeric& c);
    [[nodiscard]] Scalar arg(const ComplexNumeric& c);
    ComplexNumeric operator+(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator-(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator*(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator/(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator*(const ComplexNumeric& c, const Scalar& n);
    ComplexNumeric operator/(const ComplexNumeric& c, const Scalar& n);
    void swap(ComplexNumeric& c1, ComplexNumeric& c2) noexcept { c1.swap(c2); }
    void operator+=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 + c2; }
    void operator-=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 - c2; }
    void operator*=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 * c2; }
    void operator/=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 / c2; }
    void operator*=(ComplexNumeric& c, const Scalar& n) { c = c * n; }
    void operator/=(ComplexNumeric& c, const Scalar& n) { c = c / n; }
}

#endif