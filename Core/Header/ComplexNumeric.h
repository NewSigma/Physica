/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_COMPLEXNUMERIC_H
#define PHYSICA_COMPLEXNUMERIC_H

#include "Numerical.h"

namespace Physica::Core {
    class ComplexNumeric : private Numerical {
        Numerical imagine;
    public:
        ComplexNumeric() = default;
        ComplexNumeric(Numerical n1, Numerical n2) noexcept;
        ComplexNumeric(const ComplexNumeric& c) = default;
        ComplexNumeric(ComplexNumeric&& c) noexcept;
        /* Operators */
        friend std::ostream& operator<<(std::ostream& os, const ComplexNumeric& n);
        ComplexNumeric& operator=(const ComplexNumeric& c);
        ComplexNumeric& operator=(ComplexNumeric&& c) noexcept;
        /* Helpers */
        void swap(ComplexNumeric& c) noexcept;
        /* Getters */
        [[nodiscard]] const Numerical& getReal() const { return *this; }
        [[nodiscard]] const Numerical& getImagine() const { return imagine; }
    };
    [[nodiscard]] Numerical norm(const ComplexNumeric& c);
    [[nodiscard]] Numerical arg(const ComplexNumeric& c);
    ComplexNumeric operator+(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator-(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator*(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator/(const ComplexNumeric& c1, const ComplexNumeric& c2);
    ComplexNumeric operator*(const ComplexNumeric& c, const Numerical& n);
    ComplexNumeric operator/(const ComplexNumeric& c, const Numerical& n);
    void swap(ComplexNumeric& c1, ComplexNumeric& c2) noexcept { c1.swap(c2); }
    void operator+=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 + c2; }
    void operator-=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 - c2; }
    void operator*=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 * c2; }
    void operator/=(ComplexNumeric& c1, const ComplexNumeric& c2) { c1 = c1 / c2; }
    void operator*=(ComplexNumeric& c, const Numerical& n) { c = c * n; }
    void operator/=(ComplexNumeric& c, const Numerical& n) { c = c / n; }
}

#endif