/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include <iosfwd>
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Utils/Container/CStyleArray/CStyleArray.h"

#include "VectorExpression.h"

namespace Physica::Core {
    /*!
     * T must be either Scalar or ComplexScalar.
     */
    template<class T, size_t maxLength>
    class Vector : public Utils::CStyleArray<T, maxLength> {
        typedef Utils::CStyleArray<T, maxLength> Base;
        typedef T ScalarType;
    public:
        Vector();
        template<VectorExpressionType type, class T1, class T2>
        Vector(const VectorExpression<type, T1, T2>& expression); //NOLINT Implicit conversions is permitted.
        explicit Vector(size_t length);
        explicit Vector(const Base& array);
        explicit Vector(Base&& array) noexcept;
        Vector(std::initializer_list<T> list);
        Vector(const Vector<T, maxLength>& vec);
        Vector(Vector<T, maxLength>&& vec) noexcept;
        ~Vector() = default;
        /* Operators */
        Vector<T, maxLength>& operator=(const Vector<T, maxLength>& v) noexcept { Base::operator=(v); return *this; }
        Vector<T, maxLength>& operator=(Vector<T, maxLength>&& v) noexcept { Base::operator=(std::move(v)); return *this; }
        template<VectorExpressionType type, class T1, class T2>
        Vector<T, maxLength>& operator=(const VectorExpression<type, T1, T2>& exp);
        /* Vector Operations */
        Vector& toOpposite();
        [[nodiscard]] T toNorm() const;
        void toUnit();
        /* Getters */
        [[nodiscard]] bool isZero() const;
        /* Helpers */
        static Vector<T, Utils::Dynamic> zeroVector(size_t len);
        static Vector<T, Utils::Dynamic> randomVector(size_t len);
        static Vector simplyMultiply(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2);
    };
    /* Operators */
    template<class T, size_t maxLength>
    std::ostream& operator<<(std::ostream& os, const Vector<T, maxLength>& v);

    template<class T, size_t maxLength>
    inline VectorExpression<VectorExpressionType::Minus, Vector<T, maxLength>> operator-(const Vector<T, maxLength>& v);

    template<class T, size_t maxLength>
    inline VectorExpression<VectorExpressionType::Add, Vector<T, maxLength>, Vector<T, maxLength>>
            operator+(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2);

    template<class T, size_t maxLength>
    inline VectorExpression<VectorExpressionType::Sub, Vector<T, maxLength>, Vector<T, maxLength>>
            operator-(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2);

    template<class T, size_t maxLength>
    T operator*(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2);

    template<class T, size_t maxLength>
    VectorExpression<VectorExpressionType::Add, Vector<T, maxLength>, T> operator+(const Vector<T, maxLength>& v, const T& s);

    template<class T, size_t maxLength>
    VectorExpression<VectorExpressionType::Sub, Vector<T, maxLength>, T> operator-(const Vector<T, maxLength>& v, const T& s);

    template<class T, size_t maxLength>
    VectorExpression<VectorExpressionType::Mul, Vector<T, maxLength>, T> operator*(const Vector<T, maxLength>& v, const T& s);

    template<class T, size_t maxLength>
    VectorExpression<VectorExpressionType::Div, Vector<T, maxLength>, T> operator/(const Vector<T, maxLength>& v, const T& s);
    /* Inline Implements */
    template<class T, size_t maxLength>
    inline void operator+=(Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2) { v1 = v1 + v2; }

    template<class T, size_t maxLength>
    inline void operator-=(Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2) { v1 = v1 - v2; }

    template<class T, size_t maxLength>
    inline void operator+=(Vector<T, maxLength>& v, const T& n) { v = v + n; }

    template<class T, size_t maxLength>
    inline void operator-=(Vector<T, maxLength>& v, const T& n) { v = v - n; }

    template<class T, size_t maxLength>
    inline void operator*=(Vector<T, maxLength>& v, const T& n) { v = v * n; }

    template<class T, size_t maxLength>
    inline void operator/=(Vector<T, maxLength>& v, const T& n) { v = v / n; }

    template<class T, size_t maxLength, VectorExpressionType type, class T1, class T2>
    void operator+=(Vector<T, maxLength>& v, VectorExpression<type, T1, T2> expression) {
        const size_t length = expression.getLength();
        for(size_t i = 0; i < length; ++i)
            v[i] += expression[i];
    }

    template<class T, size_t maxLength, VectorExpressionType type, class T1, class T2>
    void operator-=(Vector<T, maxLength>& v, VectorExpression<type, T1, T2> expression) {
        const size_t length = expression.getLength();
        for(size_t i = 0; i < length; ++i)
            v[i] += expression[i];
    }

    template<class T, size_t maxLength, VectorExpressionType type, class T1, class T2>
    void operator*=(Vector<T, maxLength>& v, VectorExpression<type, T1, T2> expression) {
        const size_t length = expression.getLength();
        for(size_t i = 0; i < length; ++i)
            v[i] *= expression[i];
    }

    template<class T, size_t maxLength, VectorExpressionType type, class T1, class T2>
    void operator/=(Vector<T, maxLength>& v, VectorExpression<type, T1, T2> expression) {
        const size_t length = expression.getLength();
        for(size_t i = 0; i < length; ++i)
            v[i] /= expression[i];
    }
}

#include "VectorImpl.h"

#endif
