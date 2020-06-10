/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include <iomanip>
#include <utility>
#include "Core/Header/ComplexNumeric.h"

namespace Physica::Core {
    ComplexNumeric::ComplexNumeric(Numerical n1, Numerical n2) noexcept
    : Numerical(std::move(n1)), imagine(std::move(n2)) {}

    ComplexNumeric::ComplexNumeric(ComplexNumeric&& c) noexcept
    : Numerical(std::move(c)), imagine(std::move(c.imagine)) {}

    ComplexNumeric& ComplexNumeric::operator=(const ComplexNumeric& c) {
        if(this == &c)
            return *this;
        Numerical::operator=(c);
        imagine = c.imagine;
        return *this;
    }

    ComplexNumeric& ComplexNumeric::operator=(ComplexNumeric&& c) noexcept {
        Numerical::operator=(static_cast<Numerical&&>(c));
        imagine = std::move(c.imagine);
        return *this;
    }

    std::ostream& operator<<(std::ostream& os, const ComplexNumeric& n) {
        return os << std::setprecision(10) << double(n.getReal()) << " + i" << double(n.imagine) << std::setprecision(6);
    }

    void ComplexNumeric::swap(ComplexNumeric& c) noexcept {
        Numerical::swap(c);
        Physica::Core::swap(imagine, c.imagine);
    }

    Numerical norm(const ComplexNumeric& c) {
        return sqrt(square(c.getReal()) + square(c.getImagine()));
    }

    Numerical arg(const ComplexNumeric& c) {
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

    ComplexNumeric operator*(const ComplexNumeric& c, const Numerical& n) {
        return ComplexNumeric(c.getReal() * n, c.getImagine() * n);
    }

    ComplexNumeric operator/(const ComplexNumeric& c, const Numerical& n) {
        return ComplexNumeric(c.getReal() / n, c.getImagine() / n);
    }
}