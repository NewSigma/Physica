/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_VECTORIMPL_H
#define PHYSICA_VECTORIMPL_H
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength>::Vector() : Base() {}

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength>::Vector(const Base& array) : Base(array) {}

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength>::Vector(Base&& array) noexcept : Base(std::move(array)) {}

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength>::Vector(Vector&& vec) noexcept : Base(std::move(vec))  {}

    template<class T, size_t length, size_t maxLength>
    T Vector<T, length, maxLength>::toNorm() const {
        auto norm = T::getZero();
        for(size_t i = 0; i < CStyleArray<T, length, maxLength>::getLength(); ++i)
            norm += square(CStyleArray<T, length, maxLength>::operator[](i));
        return sqrt(norm);
    }

    template<class T, size_t length, size_t maxLength>
    void Vector<T, length, maxLength>::toUnit() {
        if(isZero())
            return;
        T norm = toNorm();
        for(size_t i = 0; i < CStyleArray<T, length, maxLength>::getLength(); ++i)
            CStyleArray<T, length, maxLength>::operator[](i) /= norm;
    }

    template<class T, size_t length, size_t maxLength>
    bool Vector<T, length, maxLength>::isZero() const {
        const auto len = CStyleArray<T, length, maxLength>::getLength();
        if(len == 0)
            return false;
        for(size_t i = 0; i < len; ++i) {
            if(!(*this)[i].isZero())
                return false;
        }
        return true;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, Dynamic, Dynamic> Vector<T, length, maxLength>::zeroVector(size_t len) {
        Vector<T, Dynamic, Dynamic> result(CStyleArray<T, Dynamic, Dynamic>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(T::getZero(), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, Dynamic, Dynamic> Vector<T, length, maxLength>::randomVector(size_t len) {
        Vector<T, Dynamic, Dynamic> result(CStyleArray<T, Dynamic, Dynamic>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(randomScalar<T::getType(), T::getErrorTrack()>(), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> Vector<T, length, maxLength>::simplyMultiply(
            const Vector<T, length, maxLength>& v1, const Vector<T, length, maxLength>& v2) {
        const auto len = v1.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v1[i] * v2[i], i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    std::ostream& operator<<(std::ostream& os, const Vector<T, length, maxLength>& v) {
        os << '(';
        const auto length_1 = v.getLength() - 1;
        for(size_t i = 0; i < length_1; ++i)
            os << v[i] << ", ";
        os << v[length_1] << ')';
        return os;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> operator-(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(-v[i], i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> operator+(const Vector<T, length, maxLength>& v1, const Vector<T, length, maxLength>& v2) {
        const auto len = v1.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v1[i] + v2[i], i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> operator-(const Vector<T, length, maxLength>& v1, const Vector<T, length, maxLength>& v2) {
        const auto len = v1.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v1[i] - v2[i], i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    T operator*(const Vector<T, length, maxLength>& v1, const Vector<T, length, maxLength>& v2) {
        const auto len = v1.getLength();
        auto result = T::getZero();
        for(size_t i = 0; i < len; ++i)
            result += v1[i] * v2[i];
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> operator+(const Vector<T, length, maxLength>& v, const T& n) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v[i] + n, i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> operator-(const Vector<T, length, maxLength>& v, const T& n) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v[i] - n, i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> operator*(const Vector<T, length, maxLength>& v, const T& n) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v[i] * n, i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> operator/(const Vector<T, length, maxLength>& v, const T& n) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v[i] / n, i);
        return result;
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> reciprocal(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(reciprocal(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> sqrt(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(sqrt(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> factorial(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(factorial(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> ln(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(ln(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> log(const Vector<T, length, maxLength>& v, const T& a) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(log(v[i], a), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> exp(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(exp(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> pow(const Vector<T, length, maxLength>& v, const T& a) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v[i] ^ a, i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> cos(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(cos(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> sin(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(sin(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> tan(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(tan(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> sec(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(sec(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> csc(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(csc(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> cot(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(cot(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arccos(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccos(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arcsin(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arcsin(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arctan(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arctan(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arcsec(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arcsec(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arccsc(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccsc(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arccot(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccot(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> cosh(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(cosh(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> sinh(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(sinh(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> tanh(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(tanh(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> sech(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(sech(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> csch(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(csch(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> coth(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(coth(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arccosh(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccosh(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arcsinh(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arcsinh(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arctanh(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arctanh(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arcsech(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arcsech(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arccsch(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccsch(v[i]), i);
        return result;
    }

    template<class T, size_t length, size_t maxLength>
    Vector<T, length, maxLength> arccoth(const Vector<T, length, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, length, maxLength> result(CStyleArray<T, length, maxLength>(len, len));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccoth(v[i]), i);
        return result;
    }
}

#endif