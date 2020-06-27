/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include <iomanip>
#include <utility>
#include "Physica/Core/Math/Algebra/ComplexNumeric.h"
#include "Physica/Core/MultiPrecition/ElementaryFunction.h"

namespace Physica::Core {
    ComplexNumeric::ComplexNumeric(MultiScalar n1, MultiScalar n2) noexcept
    : MultiScalar(std::move(n1)), imagine(std::move(n2)) {}

    ComplexNumeric::ComplexNumeric(ComplexNumeric&& c) noexcept
    : MultiScalar(std::move(c)), imagine(std::move(c.imagine)) {}

    ComplexNumeric& ComplexNumeric::operator=(const ComplexNumeric& c) {
        if(this == &c)
            return *this;
        MultiScalar::operator=(c);
        imagine = c.imagine;
        return *this;
    }

    ComplexNumeric& ComplexNumeric::operator=(ComplexNumeric&& c) noexcept {
        MultiScalar::operator=(static_cast<MultiScalar&&>(c));
        imagine = std::move(c.imagine);
        return *this;
    }

    std::ostream& operator<<(std::ostream& os, const ComplexNumeric& n) {
        return os << std::setprecision(10) << double(n.getReal()) << " + i" << double(n.imagine) << std::setprecision(6);
    }

    void ComplexNumeric::swap(ComplexNumeric& c) noexcept {
        MultiScalar::swap(c);
        Physica::Core::swap(imagine, c.imagine);
    }

    MultiScalar norm(const ComplexNumeric& c) {
        return sqrt(square(c.getReal()) + square(c.getImagine()));
    }

    MultiScalar arg(const ComplexNumeric& c) {
        return arctan(c.getImagine() / c.getReal());
    }

    ComplexNumeric operator+(const ComplexNumeric& c1, const ComplexNumeric& c2) {
        return ComplexNumeric(c1.getReal() + c2.getReal(), c1.getImagine() + c2.getImagine());
    }

    ComplexNumeric operator-(const ComplexNumeric& c1, const ComplexNumeric& c2) {
        return ComplexNumeric(c1.getReal() - c2.getReal(), c1.getImagine() - c2.getImagine());
    }

    ComplexNumeric operator*(const ComplexNumeric& c1, const ComplexNumeric& c2) {
        return ComplexNumeric(c1.getReal() * c2.getReal() - c1.getImagine() * c2.getImagine()
                , c1.getReal() * c2.getImagine() + c1.getImagine() * c2.getReal());
    }

    ComplexNumeric operator/(const ComplexNumeric& c1, const ComplexNumeric& c2) {
        return ComplexNumeric(c1.getReal() * c2.getReal() + c1.getImagine() * c2.getImagine()
                ,  c1.getImagine() * c2.getReal() - c1.getReal() * c2.getImagine())
                / (square(c2.getReal()) + square(c2.getImagine()));
    }

    ComplexNumeric operator*(const ComplexNumeric& c, const MultiScalar& n) {
        return ComplexNumeric(c.getReal() * n, c.getImagine() * n);
    }

    ComplexNumeric operator/(const ComplexNumeric& c, const MultiScalar& n) {
        return ComplexNumeric(c.getReal() / n, c.getImagine() / n);
    }
}