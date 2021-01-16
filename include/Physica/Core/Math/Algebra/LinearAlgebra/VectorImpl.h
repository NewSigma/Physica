/*
 * Copyright 2020-2021 WeiBo He.
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
#pragma once

#include <cassert>
/**
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    template<class T, size_t Length, size_t maxLength>
    template<VectorExpressionType type, class T1, class T2>
    Vector<T, Length, maxLength>::Vector(const VectorExpression<type, T1, T2>& expression) : Base(expression.getLength()) {
        const size_t length = expression.getLength();
        for (size_t i = 0; i < length; ++i)
            Base::init(expression[i], i);
        Base::setLength(length);
    }

    template<class T, size_t Length, size_t maxLength>
    template<VectorExpressionType type, class T1, class T2>
    Vector<T, Length, maxLength>& Vector<T, Length, maxLength>::operator=(const VectorExpression<type, T1, T2>& exp) {
        Base::clear();
        const size_t length = exp.getLength();
        for (size_t i = 0; i < length; ++i)
            Base::init(exp[i], i);
        Base::setLength(length);
        return *this;
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength>& Vector<T, Length, maxLength>::toOpposite() {
        const auto end = Base::end();
        for (auto ite = Base::begin(); ite != end; ++ite)
            (*ite).toOpposite();
        return *this;
    }

    template<class T, size_t Length, size_t maxLength>
    T Vector<T, Length, maxLength>::toNorm() const {
        auto norm = T::getZero();
        for(auto ite = Base::cbegin(); ite != Base::cend(); ++ite)
            norm += square(*ite);
        return sqrt(norm);
    }

    template<class T, size_t Length, size_t maxLength>
    void Vector<T, Length, maxLength>::toUnit() {
        T norm = toNorm();
        if (norm.isZero())
            return;
        const auto end = Base::end();
        for (auto ite = Base::begin(); ite != end; ++ite)
            (*ite) /= norm;
    }

    template<class T, size_t Length, size_t maxLength>
    bool Vector<T, Length, maxLength>::isZero() const {
        const auto length = Base::getLength();
        if(length == 0)
            return false;
        for(size_t i = 0; i < length; ++i) {
            if(!(*this)[i].isZero())
                return false;
        }
        return true;
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T> Vector<T, Length, maxLength>::zeroVector(size_t len) {
        Vector<T> result(len);
        for(size_t i = 0; i < len; ++i)
            result.allocate(T::getZero(), i);
        result.setLength(len);
        return result;
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T> Vector<T, Length, maxLength>::randomVector(size_t len) {
        Vector<T> result(len);
        for (size_t i = 0; i < len; ++i)
            result.allocate(randomScalar<T::getType(), T::getErrorTrack()>(), i);
        result.setLength(len);
        return result;
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> Vector<T, Length, maxLength>::simplyMultiply(const Vector& v1, const Vector& v2) {
        const auto len = v1.getLength();
        assert(len == v2.getLength());
        Vector<T, Length, maxLength> result(len);
        for (size_t i = 0; i < len; ++i)
            result.init(v1[i] * v2[i], i);
        result.setLength(len);
        return result;
    }

    template<class T, size_t Length, size_t maxLength>
    std::ostream& operator<<(std::ostream& os, const Vector<T, Length, maxLength>& v) {
        os << '(';
        size_t length = v.getLength();
        if (length) {
            --length;
            for (size_t i = 0; i < length; ++i)
                os << v[i] << ", ";
            os << v[length];
        }
        os << ')';
        return os;
    }

    template<class T, size_t Length, size_t maxLength>
    inline VectorExpression<VectorExpressionType::Minus, Vector<T, Length, maxLength>> operator-(const Vector<T, Length, maxLength>& v) {
        return VectorExpression<VectorExpressionType::Minus, Vector<T, Length, maxLength>>(v);
    }

    template<class T, size_t Length, size_t maxLength>
    inline VectorExpression<VectorExpressionType::Add, Vector<T, Length, maxLength>, Vector<T, Length, maxLength>>
            operator+(const Vector<T, Length, maxLength>& v1, const Vector<T, Length, maxLength>& v2) {
        return VectorExpression<VectorExpressionType::Add, Vector<T, Length, maxLength>, Vector<T, Length, maxLength>>(v1, v2);
    }

    template<class T, size_t Length, size_t maxLength>
    inline VectorExpression<VectorExpressionType::Sub, Vector<T, Length, maxLength>, Vector<T, Length, maxLength>>
            operator-(const Vector<T, Length, maxLength>& v1, const Vector<T, Length, maxLength>& v2) {
        return VectorExpression<VectorExpressionType::Sub, Vector<T, Length, maxLength>, Vector<T, Length, maxLength>>(v1, v2);
    }

    template<class T, size_t Length, size_t maxLength>
    T operator*(const Vector<T, Length, maxLength>& v1, const Vector<T, Length, maxLength>& v2) {
        const auto len = v1.getLength();
        assert(len == v2.getLength());
        auto result = T::getZero();
        for(size_t i = 0; i < len; ++i)
            result += v1[i] * v2[i];
        return result;
    }

    template<class T, size_t Length, size_t maxLength>
    VectorExpression<VectorExpressionType::Add, Vector<T, Length, maxLength>, T> operator+(const Vector<T, Length, maxLength>& v, const T& s) {
        return VectorExpression<VectorExpressionType::Add, Vector<T, Length, maxLength>, T>(v, s);
    }

    template<class T, size_t Length, size_t maxLength>
    VectorExpression<VectorExpressionType::Sub, Vector<T, Length, maxLength>, T> operator-(const Vector<T, Length, maxLength>& v, const T& s) {
        return VectorExpression<VectorExpressionType::Sub, Vector<T, Length, maxLength>, T>(v, s);
    }

    template<class T, size_t Length, size_t maxLength>
    VectorExpression<VectorExpressionType::Mul, Vector<T, Length, maxLength>, T> operator*(const Vector<T, Length, maxLength>& v, const T& s) {
        return VectorExpression<VectorExpressionType::Mul, Vector<T, Length, maxLength>, T>(v, s);
    }

    template<class T, size_t Length, size_t maxLength>
    VectorExpression<VectorExpressionType::Div, Vector<T, Length, maxLength>, T> operator/(const Vector<T, Length, maxLength>& v, const T& s) {
        return VectorExpression<VectorExpressionType::Div, Vector<T, Length, maxLength>, T>(v, s);
    }

    template<class T, size_t Length, size_t maxLength, VectorExpressionType type, class T1, class T2>
    void operator+=(Vector<T, Length, maxLength>& v1, const VectorExpression<type, T1, T2>& exp) {
        for (size_t i = 0; i < exp.getLength(); ++i)
            v1[i] = v1[i] + exp[i];
    }

    template<class T, size_t Length, size_t maxLength, VectorExpressionType type, class T1, class T2>
    void operator-=(Vector<T, Length, maxLength>& v1, const VectorExpression<type, T1, T2>& exp) {
        for (size_t i = 0; i < exp.getLength(); ++i)
            v1[i] = v1[i] - exp[i];
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    //Optimize: the following functions may be speed up using expression templates.
    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> reciprocal(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(reciprocal(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> sqrt(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(sqrt(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> factorial(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(factorial(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> ln(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(ln(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> log(const Vector<T, Length, maxLength>& v, const T& a) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(log(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> exp(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(exp(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> pow(const Vector<T, Length, maxLength>& v, const T& a) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(pow(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> cos(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(cos(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> sin(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(sin(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> tan(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(tan(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> sec(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(sec(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> csc(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(csc(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> cot(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(cot(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arccos(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccos(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arcsin(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arcsin(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arctan(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arctan(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arcsec(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arcsec(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arccsc(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccsc(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arccot(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccot(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> cosh(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(cosh(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> sinh(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(sinh(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> tanh(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(tanh(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> sech(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(sech(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> csch(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(csch(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> coth(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(coth(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arccosh(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccosh(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arcsinh(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arcsinh(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arctanh(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arctanh(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arcsech(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arcsech(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arccsch(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccsch(*ite), i);
    }

    template<class T, size_t Length, size_t maxLength>
    Vector<T, Length, maxLength> arccoth(const Vector<T, Length, maxLength>& v) {
        Vector<T, Length, maxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccoth(*ite), i);
    }
}
