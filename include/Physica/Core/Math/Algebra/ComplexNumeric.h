/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_COMPLEXNUMERIC_H
#define PHYSICA_COMPLEXNUMERIC_H

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    class ComplexNumeric : private MultiScalar {
        MultiScalar imagine;
    public:
        ComplexNumeric() = default;
        ComplexNumeric(MultiScalar n1, MultiScalar n2) noexcept;
        ComplexNumeric(const ComplexNumeric& c) = default;
        ComplexNumeric(ComplexNumeric&& c) noexcept;
        /* Operators */
        friend std::ostream& operator<<(std::ostream& os, const ComplexNumeric& n);
        ComplexNumeric& operator=(const ComplexNumeric& c);
        ComplexNumeric& operator=(ComplexNumeric&& c) noexcept;
        /* Helpers */
        void swap(ComplexNumeric& c) noexcept;
        /* Getters */
        [[nodiscard]] const MultiScalar& getReal() const { return *this; }
        [[nodiscard]] const MultiScalar& getImagine() const { return imagine; }
    };
    [[nodiscard]] MultiScalar norm(const ComplexNumeric& c);
    [[nodiscard]] MultiScalar arg(const ComplexNumeric& c);
    ComplexNumeric operator+(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator-(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator*(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator/(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator*(const ComplexNumeric& c, const MultiScalar& n);
    ComplexNumeric operator/(const ComplexNumeric& c, const MultiScalar& n);
    void swap(ComplexNumeric& c1, ComplexNumeric& c2) noexcept { c1.swap(c2); }
    void operator+=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 + c2; }
    void operator-=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 - c2; }
    void operator*=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 * c2; }
    void operator/=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 / c2; }
    void operator*=(ComplexNumeric& c, const MultiScalar& n) { c = c * n; }
    void operator/=(ComplexNumeric& c, const MultiScalar& n) { c = c / n; }
}

#endif