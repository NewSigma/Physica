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
    template<class Derived>
    Vector<T, Length, MaxLength>::Vector(const RValueVector<Derived>& v) : Storage(v.getLength()) {
        v.assignTo(*this);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength>& Vector<T, Length, MaxLength>::toOpposite() {
        const auto end = Storage::end();
        for (auto ite = Storage::begin(); ite != end; ++ite)
            (*ite).toOpposite();
        return *this;
    }

    template<class T, size_t Length, size_t MaxLength>
    void Vector<T, Length, MaxLength>::toUnit() {
        T norm = Base::norm();
        if (norm.isZero())
            return;
        const auto end = Storage::end();
        for (auto ite = Storage::begin(); ite != end; ++ite)
            (*ite) /= norm;
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> Vector<T, Length, MaxLength>::Zeros(size_t len) {
        Vector<T, Length, MaxLength> result{};
        result.reserve(len);
        for(size_t i = 0; i < len; ++i)
            result.get_allocator().construct(result.data() + i, T::Zero());
        result.setLength(len);
        return result;
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> Vector<T, Length, MaxLength>::randomVector(size_t len) {
        Vector<T, Length, MaxLength> result{};
        result.reserve(len);
        for (size_t i = 0; i < len; ++i)
            result.get_allocator().construct(result.data() + i, randomScalar<T>());
        result.setLength(len);
        return result;
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> Vector<T, Length, MaxLength>::randomVector(const Vector& v1, const Vector& v2) {
        assert(v1.getLength() == v2.getLength());
        Vector<T, Length, MaxLength> result = randomVector(v1.getLength());
        result = v1 + multiply((v2 - v1), result);
        return result;
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> Vector<T, Length, MaxLength>::linspace(T from, T to, size_t count) {
        assert(from < to);
        const T step = (to - from) / T(count - 1);
        Vector result = Vector(count);
        for (size_t i = 0; i < count; ++i) {
            result[i] = from;
            from += step;
        }
        return result;
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    //Optimize: the following functions may be speed up using expression templates.
    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> factorial(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(factorial(*ite), i);
    }

    template<class T, size_t Length, size_t MaxLength>
    Vector<T, Length, MaxLength> exp(const Vector<T, Length, MaxLength>& v) {
        Vector<T, Length, MaxLength> result(v.getLength());
        size_t i = 0;
        for (auto ite = v.cbegin(); ite != v.cend(); ++ite, ++i)
            result.init(exp(*ite), i);
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
