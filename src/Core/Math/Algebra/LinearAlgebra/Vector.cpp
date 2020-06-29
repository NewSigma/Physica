/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtCore/qlogging.h>
#include <cstring>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector.h"
#include "Physica/Core/MultiPrecition/ScalarImpl/ElementaryFunction.h"

namespace Physica::Core {
    Vector::Vector() : CStyleArray<MultiScalar, Dynamic>() {}

    Vector::Vector(size_t capacity) : CStyleArray<MultiScalar, Dynamic>(capacity) {}
    /*!
     * Low level api. Designed for performance.
     * Warning: The first \length elements have not allocated. DO NOT try to visit them.
     */
    Vector::Vector(size_t length, size_t capacity) : CStyleArray<MultiScalar, Dynamic>(length, capacity) {}
    //!Convenience function to create a 2D Vector.
    Vector::Vector(const MultiScalar& n1, const MultiScalar& n2) : CStyleArray<MultiScalar, Dynamic>(2, 2) {
        allocate(n1, 0);
        allocate(n2, 1);
    }
    //!Convenience function to create a 3D Vector.
    Vector::Vector(const MultiScalar& n1, const MultiScalar& n2, const MultiScalar& n3)
            : CStyleArray<MultiScalar, Dynamic>(3, 3) {
        allocate(n1, 0);
        allocate(n2, 1);
        allocate(n3, 2);
    }

    Vector::Vector(const CStyleArray<MultiScalar, Dynamic>& array) : CStyleArray<MultiScalar, Dynamic>(array) {}

    Vector::Vector(CStyleArray<MultiScalar, Dynamic>&& array) noexcept : CStyleArray<MultiScalar, Dynamic>(std::move(array)) {}

    Vector::Vector(Vector&& vec) noexcept : CStyleArray<MultiScalar, Dynamic>(std::move(vec))  {}

    std::ostream& operator<<(std::ostream& os, const Vector& v) {
        os << '(';
        const auto length_1 = v.getLength() - 1;
        for(size_t i = 0; i < length_1; ++i)
            os << double(v[i]) << ", ";
        os << double(v[length_1]) << ')';
        return os;
    }

    MultiScalar Vector::toNorm() const {
        MultiScalar norm(static_cast<SignedScalarUnit>(0));
        for(size_t i = 0; i < getLength(); ++i)
            norm += square((*this)[i]);
        return sqrt(norm);
    }

    void Vector::toUnit() {
        if(isZeroVector())
            return;
        MultiScalar norm = toNorm();
        for(size_t i = 0; i < getLength(); ++i)
            (*this)[i] /= norm;
    }

    bool Vector::isZeroVector() const {
        const auto length = getLength();
        if(length == 0)
            return false;
        for(size_t i = 0; i < length; ++i) {
            if(!(*this)[i].isZero())
                return false;
        }
        return true;
    }

    MultiScalar Vector::toArg(size_t axe) const {
        return toNorm() / (*this)[axe];
    }

    Vector Vector::zeroVector(size_t length) {
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(MultiScalar(static_cast<SignedScalarUnit>(0)));
        return result;
    }

    Vector Vector::randomVector(size_t length) {
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(randomScalar<MultiPrecision, true>());
        return result;
    }

    Vector Vector::simplyMultiply(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(v1[i] * v2[i], i);
        return result;
    }

    Vector operator-(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(-v[i], i);
        return result;
    }

    bool operator==(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        if(length != v2.getLength())
            return false;
        for(size_t i = 0; i < length; ++i)
            if(v1[i] != v2[i])
                return false;
        return true;
    }

    Vector operator+(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(v1[i] + v2[i], i);
        return result;
    }

    Vector operator-(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(v1[i] - v2[i], i);
        return result;
    }

    MultiScalar operator*(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        MultiScalar result(static_cast<SignedScalarUnit>(0));
        for(size_t i = 0; i < length; ++i)
            result += v1[i] * v2[i];
        return result;
    }

    Vector operator+(const Vector& v, const MultiScalar& n) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(v[i] + n, i);
        return result;
    }

    Vector operator-(const Vector& v, const MultiScalar& n) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(v[i] - n, i);
        return result;
    }

    Vector operator*(const Vector& v, const MultiScalar& n) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(v[i] * n, i);
        return result;
    }

    Vector operator/(const Vector& v, const MultiScalar& n) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(v[i] / n, i);
        return result;
    }
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    Vector reciprocal(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(reciprocal(v[i]), i);
        return result;
    }

    Vector sqrt(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(sqrt(v[i]), i);
        return result;
    }

    Vector factorial(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(factorial(v[i]), i);
        return result;
    }

    Vector ln(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(ln(v[i]), i);
        return result;
    }

    Vector log(const Vector& v, const MultiScalar& a) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(log(v[i], a), i);
        return result;
    }

    Vector exp(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(exp(v[i]), i);
        return result;
    }

    Vector pow(const Vector& v, const MultiScalar& a) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(v[i] ^ a, i);
        return result;
    }

    Vector cos(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(cos(v[i]), i);
        return result;
    }

    Vector sin(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(sin(v[i]), i);
        return result;
    }

    Vector tan(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(tan(v[i]), i);
        return result;
    }

    Vector sec(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(sec(v[i]), i);
        return result;
    }

    Vector csc(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(csc(v[i]), i);
        return result;
    }

    Vector cot(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(cot(v[i]), i);
        return result;
    }

    Vector arccos(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccos(v[i]), i);
        return result;
    }

    Vector arcsin(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arcsin(v[i]), i);
        return result;
    }

    Vector arctan(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arctan(v[i]), i);
        return result;
    }

    Vector arcsec(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arcsec(v[i]), i);
        return result;
    }

    Vector arccsc(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccsc(v[i]), i);
        return result;
    }

    Vector arccot(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccot(v[i]), i);
        return result;
    }

    Vector cosh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(cosh(v[i]), i);
        return result;
    }

    Vector sinh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(sinh(v[i]), i);
        return result;
    }

    Vector tanh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(tanh(v[i]), i);
        return result;
    }

    Vector sech(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(sech(v[i]), i);
        return result;
    }

    Vector csch(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(csch(v[i]), i);
        return result;
    }

    Vector coth(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(coth(v[i]), i);
        return result;
    }

    Vector arccosh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccosh(v[i]), i);
        return result;
    }

    Vector arcsinh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arcsinh(v[i]), i);
        return result;
    }

    Vector arctanh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arctanh(v[i]), i);
        return result;
    }

    Vector arcsech(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arcsech(v[i]), i);
        return result;
    }

    Vector arccsch(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccsch(v[i]), i);
        return result;
    }

    Vector arccoth(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length, length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccoth(v[i]), i);
        return result;
    }
}