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
    template<class T, size_t maxLength>
    Vector<T, maxLength>::Vector() : Base() {}
    /*!
     * Elements must not be accessed before they are initialized.
     */
    template<class T, size_t maxLength>
    Vector<T, maxLength>::Vector(size_t length) : Base(length) {}

    template<class T, size_t maxLength>
    Vector<T, maxLength>::Vector(const Base& array) : Base(array) {}

    template<class T, size_t maxLength>
    Vector<T, maxLength>::Vector(Base&& array) noexcept : Base(std::move(array)) {}

    template<class T, size_t maxLength>
    Vector<T, maxLength>::Vector(const Vector<T, maxLength>& vec) : Base(vec) {}

    template<class T, size_t maxLength>
    Vector<T, maxLength>::Vector(Vector&& vec) noexcept : Base(std::move(vec))  {}

    template<class T, size_t maxLength>
    Vector<T, maxLength>& Vector<T, maxLength>::toOpposite() {
        const auto length = Base::getLength();
        for(size_t i = 0; i < length; ++i)
            Base::operator[](i).toOpposite();
    }

    template<class T, size_t maxLength>
    T Vector<T, maxLength>::toNorm() const {
        auto norm = T::getZero();
        for(size_t i = 0; i < CStyleArray<T, maxLength>::getLength(); ++i)
            norm += square(CStyleArray<T, maxLength>::operator[](i));
        return sqrt(norm);
    }

    template<class T, size_t maxLength>
    void Vector<T, maxLength>::toUnit() {
        if(isZero())
            return;
        T norm = toNorm();
        for(size_t i = 0; i < CStyleArray<T, maxLength>::getLength(); ++i)
            CStyleArray<T, maxLength>::operator[](i) /= norm;
    }

    template<class T, size_t maxLength>
    bool Vector<T, maxLength>::isZero() const {
        const auto len = CStyleArray<T, maxLength>::getLength();
        if(len == 0)
            return false;
        for(size_t i = 0; i < len; ++i) {
            if(!(*this)[i].isZero())
                return false;
        }
        return true;
    }

    template<class T, size_t maxLength>
    Vector<T, Dynamic> Vector<T, maxLength>::zeroVector(size_t len) {
        Vector<T, Dynamic> result((CStyleArray<T, Dynamic>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(T::getZero(), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, Dynamic> Vector<T, maxLength>::randomVector(size_t len) {
        Vector<T, Dynamic> result((CStyleArray<T, Dynamic>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(randomScalar<T::getType(), T::getErrorTrack()>(), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> Vector<T, maxLength>::simplyMultiply(
            const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2) {
        const auto len = v1.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v1[i] * v2[i], i);
        return result;
    }

    template<class T, size_t maxLength>
    std::ostream& operator<<(std::ostream& os, const Vector<T, maxLength>& v) {
        os << '(';
        const auto length_1 = v.getLength() - 1;
        for(size_t i = 0; i < length_1; ++i)
            os << v[i] << ", ";
        os << v[length_1] << ')';
        return os;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator-(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(-v[i], i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator+(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2) {
        const auto len = v1.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v1[i] + v2[i], i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator-(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2) {
        const auto len = v1.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v1[i] - v2[i], i);
        return result;
    }

    template<class T, size_t maxLength>
    T operator*(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2) {
        const auto len = v1.getLength();
        auto result = T::getZero();
        for(size_t i = 0; i < len; ++i)
            result += v1[i] * v2[i];
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator+(const Vector<T, maxLength>& v, const T& n) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v[i] + n, i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator-(const Vector<T, maxLength>& v, const T& n) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v[i] - n, i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator*(const Vector<T, maxLength>& v, const T& n) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v[i] * n, i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator/(const Vector<T, maxLength>& v, const T& n) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v[i] / n, i);
        return result;
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class T, size_t maxLength>
    Vector<T, maxLength> reciprocal(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(reciprocal(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> sqrt(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(sqrt(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> factorial(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(factorial(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> ln(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(ln(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> log(const Vector<T, maxLength>& v, const T& a) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(log(v[i], a), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> exp(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(exp(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> pow(const Vector<T, maxLength>& v, const T& a) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(v[i] ^ a, i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> cos(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(cos(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> sin(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(sin(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> tan(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(tan(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> sec(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(sec(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> csc(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(csc(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> cot(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(cot(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arccos(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccos(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arcsin(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arcsin(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arctan(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arctan(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arcsec(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arcsec(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arccsc(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccsc(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arccot(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccot(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> cosh(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(cosh(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> sinh(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(sinh(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> tanh(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(tanh(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> sech(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(sech(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> csch(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(csch(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> coth(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(coth(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arccosh(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccosh(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arcsinh(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arcsinh(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arctanh(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arctanh(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arcsech(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arcsech(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arccsch(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccsch(v[i]), i);
        return result;
    }

    template<class T, size_t maxLength>
    Vector<T, maxLength> arccoth(const Vector<T, maxLength>& v) {
        const auto len = v.getLength();
        Vector<T, maxLength> result((CStyleArray<T, maxLength>(len)));
        for(size_t i = 0; i < len; ++i)
            result.allocate(arccoth(v[i]), i);
        return result;
    }
}

#endif