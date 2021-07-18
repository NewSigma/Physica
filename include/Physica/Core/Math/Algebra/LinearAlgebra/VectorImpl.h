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
    template<class T, size_t Length, size_t MaxLength>
    template<VectorExpressionType type, class T1, class T2>
    Vector<T, Length, MaxLength>::Vector(const VectorExpression<type, T1, T2>& expression) : Base(expression.getLength()) {
        const size_t length = expression.getLength();
        for (size_t i = 0; i < length; ++i)
            Base::init(expression[i], i);
        Base::setLength(length);
    }

    template<class T, size_t Length, size_t MaxLength>
    template<VectorExpressionType type, class T1, class T2>
    Vector<T, Length, MaxLength>& Vector<T, Length, MaxLength>::operator=(const VectorExpression<type, T1, T2>& exp) {
        Base::operator=(exp.calc());
        return *this;
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength>& Vector<T, Length, MaxLength>::toOpposite() {
        const auto end = Base::end();
        for (auto ite = Base::begin(); ite != end; ++ite)
            (*ite).toOpposite();
        return *this;
    }

    template<class T, size_t Length, size_t MaxLength>
    void Vector<T, Length, MaxLength>::toUnit() {
        T norm = norm();
        if (norm.isZero())
            return;
        const auto end = Base::end();
        for (auto ite = Base::begin(); ite != end; ++ite)
            (*ite) /= norm;
    }

    template<class T, size_t Length, size_t MaxLength>
    typename Vector<T, Length, MaxLength>::ColMatrix Vector<T, Length, MaxLength>::copyToColMatrix() const {
        ColMatrix mat(this->getLength(), 1);
        mat[0] = *this;
        return mat;
    }

    template<class T, size_t Length, size_t MaxLength>
    typename Vector<T, Length, MaxLength>::ColMatrix Vector<T, Length, MaxLength>::moveToColMatrix() {
        ColMatrix mat(this->getLength(), 1);
        mat[0] = std::move(*this);
        return mat;
    }

    template<class T, size_t Length, size_t MaxLength>
    typename Vector<T, Length, MaxLength>::RowMatrix Vector<T, Length, MaxLength>::copyToRowMatrix() const {
        RowMatrix mat(1, this->getLength());
        mat[0] = std::move(*this);
        return mat;
    }

    template<class T, size_t Length, size_t MaxLength>
    typename Vector<T, Length, MaxLength>::RowMatrix Vector<T, Length, MaxLength>::moveToRowMatrix() {
        RowMatrix mat(1, this->getLength());
        mat[0] = std::move(*this);
        return mat;
    }

    template<class T, size_t Length, size_t MaxLength>
    bool Vector<T, Length, MaxLength>::isZero() const {
        const auto length = Base::getLength();
        if(length == 0)
            return false;
        for(size_t i = 0; i < length; ++i) {
            if(!(*this)[i].isZero())
                return false;
        }
        return true;
    }

    template<class T, size_t Length, size_t MaxLength>
    T Vector<T, Length, MaxLength>::max() const {
        assert(Base::getLength() != 0);
        T result((*this)[0]);
        for (size_t i = 1; i < Base::getLength(); ++i) {
            if ((*this)[i] > result)
                result = (*this)[i];
        }
        return result;
    }

    template<class T, size_t Length, size_t MaxLength>
    T Vector<T, Length, MaxLength>::min() const {
        assert(Base::getLength() != 0);
        T result((*this)[0]);
        for (size_t i = 1; i < Base::getLength(); ++i) {
            if ((*this)[i] < result)
                result = (*this)[i];
        }
        return result;
    }

    template<class T, size_t Length, size_t MaxLength>
    T Vector<T, Length, MaxLength>::norm() const {
        return sqrt(squaredNorm());
    }

    template<class T, size_t Length, size_t MaxLength>
    T Vector<T, Length, MaxLength>::squaredNorm() const {
        auto result = T::getZero();
        for(auto ite = Base::cbegin(); ite != Base::cend(); ++ite)
            result += square(*ite);
        return result;
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T> Vector<T, Length, MaxLength>::zeroVector(size_t len) {
        Vector<T> result(len);
        for(size_t i = 0; i < len; ++i)
            result.allocate(T::getZero(), i);
        result.setLength(len);
        return result;
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T> Vector<T, Length, MaxLength>::randomVector(size_t len) {
        Vector<T> result(len);
        for (size_t i = 0; i < len; ++i)
            result.allocate(randomScalar<T::getType(), T::getErrorTrack()>(), i);
        result.setLength(len);
        return result;
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> Vector<T, Length, MaxLength>::simplyMultiply(const Vector& v1, const Vector& v2) {
        const auto len = v1.getLength();
        assert(len == v2.getLength());
        Vector<T, Length, MaxLength> result(len);
        for (size_t i = 0; i < len; ++i)
            result.init(v1[i] * v2[i], i);
        result.setLength(len);
        return result;
    }

    template<class T, size_t Length, size_t MaxLength>
    std::ostream& operator<<(std::ostream& os, const Vector<T, Length, MaxLength>& v) {
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

    template<class T, size_t Length, size_t MaxLength>
    inline VectorExpression<VectorExpressionType::Minus, Vector<T, Length, MaxLength>> operator-(const Vector<T, Length, MaxLength>& v) {
        return VectorExpression<VectorExpressionType::Minus, Vector<T, Length, MaxLength>>(v);
    }

    template<class T, size_t Length, size_t MaxLength>
    inline VectorExpression<VectorExpressionType::Add, Vector<T, Length, MaxLength>, Vector<T, Length, MaxLength>>
            operator+(const Vector<T, Length, MaxLength>& v1, const Vector<T, Length, MaxLength>& v2) {
        return VectorExpression<VectorExpressionType::Add, Vector<T, Length, MaxLength>, Vector<T, Length, MaxLength>>(v1, v2);
    }

    template<class T, size_t Length, size_t MaxLength>
    inline VectorExpression<VectorExpressionType::Sub, Vector<T, Length, MaxLength>, Vector<T, Length, MaxLength>>
            operator-(const Vector<T, Length, MaxLength>& v1, const Vector<T, Length, MaxLength>& v2) {
        return VectorExpression<VectorExpressionType::Sub, Vector<T, Length, MaxLength>, Vector<T, Length, MaxLength>>(v1, v2);
    }

    template<class T, size_t Length, size_t MaxLength>
    T operator*(const Vector<T, Length, MaxLength>& v1, const Vector<T, Length, MaxLength>& v2) {
        const auto len = v1.getLength();
        assert(len == v2.getLength());
        auto result = T::getZero();
        for(size_t i = 0; i < len; ++i)
            result += v1[i] * v2[i];
        return result;
    }

    template<class T, size_t Length, size_t MaxLength>
    VectorExpression<VectorExpressionType::Add, Vector<T, Length, MaxLength>, T> operator+(const Vector<T, Length, MaxLength>& v, const T& s) {
        return VectorExpression<VectorExpressionType::Add, Vector<T, Length, MaxLength>, T>(v, s);
    }

    template<class T, size_t Length, size_t MaxLength>
    VectorExpression<VectorExpressionType::Sub, Vector<T, Length, MaxLength>, T> operator-(const Vector<T, Length, MaxLength>& v, const T& s) {
        return VectorExpression<VectorExpressionType::Sub, Vector<T, Length, MaxLength>, T>(v, s);
    }

    template<class T, size_t Length, size_t MaxLength>
    VectorExpression<VectorExpressionType::Mul, Vector<T, Length, MaxLength>, T> operator*(const Vector<T, Length, MaxLength>& v, const T& s) {
        return VectorExpression<VectorExpressionType::Mul, Vector<T, Length, MaxLength>, T>(v, s);
    }

    template<class T, size_t Length, size_t MaxLength>
    VectorExpression<VectorExpressionType::Div, Vector<T, Length, MaxLength>, T> operator/(const Vector<T, Length, MaxLength>& v, const T& s) {
        return VectorExpression<VectorExpressionType::Div, Vector<T, Length, MaxLength>, T>(v, s);
    }

    template<class T, size_t Length, size_t MaxLength, VectorExpressionType type, class T1, class T2>
    void operator+=(Vector<T, Length, MaxLength>& v1, const VectorExpression<type, T1, T2>& exp) {
        for (size_t i = 0; i < exp.getLength(); ++i)
            v1[i] = v1[i] + exp[i];
    }

    template<class T, size_t Length, size_t MaxLength, VectorExpressionType type, class T1, class T2>
    void operator-=(Vector<T, Length, MaxLength>& v1, const VectorExpression<type, T1, T2>& exp) {
        for (size_t i = 0; i < exp.getLength(); ++i)
            v1[i] = v1[i] - exp[i];
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    //Optimize: the following functions may be speed up using expression templates.
    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> reciprocal(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(reciprocal(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> sqrt(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(sqrt(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> factorial(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(factorial(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> ln(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(ln(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> log(const Vector<T, Length, MaxLength>& v, const T& a) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(log(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> exp(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(exp(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> pow(const Vector<T, Length, MaxLength>& v, const T& a) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(pow(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> cos(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(cos(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> sin(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(sin(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> tan(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(tan(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> sec(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(sec(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> csc(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(csc(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> cot(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(cot(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arccos(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccos(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arcsin(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arcsin(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arctan(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arctan(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arcsec(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arcsec(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arccsc(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccsc(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arccot(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccot(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> cosh(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(cosh(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> sinh(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(sinh(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> tanh(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(tanh(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> sech(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(sech(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> csch(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(csch(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> coth(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(coth(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arccosh(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccosh(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arcsinh(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arcsinh(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arctanh(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arctanh(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arcsech(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arcsech(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arccsch(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccsch(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> arccoth(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(arccoth(*ite), i);
    }
}
