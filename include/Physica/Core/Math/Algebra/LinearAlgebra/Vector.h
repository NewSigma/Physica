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

#include <iosfwd>
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Utils/Container/Array/Array.h"
#include "VectorExpression.h"

namespace Physica::Core {
    /**
     * T must be either Scalar or ComplexScalar.
     * 
     * Default template arguments are defined in \file VectorExpression.h
     */
    template<class T, size_t Length, size_t maxLength>
    class Vector : public Utils::Array<T, Length, maxLength> {
    public:
        using ScalarType = T;
    private:
        static_assert(Length == Dynamic || Length == maxLength, "maxLength of fixed vector must equals to its length.");
        using Base = Utils::Array<T, Length, maxLength>;
    public:
        using Base::Base;
        Vector() = default;
        template<VectorExpressionType type, class T1, class T2>
        Vector(const VectorExpression<type, T1, T2>& expression); //NOLINT Implicit conversions is permitted.
        Vector(const Vector&) = default;
        Vector(Vector&&) noexcept = default;
        ~Vector() = default;
        /* Operators */
        Vector& operator=(const Vector&) = default;
        Vector& operator=(Vector&&) noexcept = default;
        template<VectorExpressionType type, class T1, class T2>
        Vector& operator=(const VectorExpression<type, T1, T2>& exp);
        /* Vector Operations */
        Vector& toOpposite();
        [[nodiscard]] T toNorm() const;
        void toUnit();
        /* Getters */
        [[nodiscard]] bool isZero() const;
        /* Helpers */
        static Vector<T> zeroVector(size_t len);
        static Vector<T> randomVector(size_t len);
        static Vector simplyMultiply(const Vector& v1, const Vector& v2);
    };
    /* Operators */
    template<class T, size_t Length, size_t maxLength>
    std::ostream& operator<<(std::ostream& os, const Vector<T, Length, maxLength>& v);

    template<class T, size_t Length, size_t maxLength>
    inline VectorExpression<VectorExpressionType::Minus, Vector<T, Length, maxLength>> operator-(const Vector<T, Length, maxLength>& v);

    template<class T, size_t Length, size_t maxLength>
    inline VectorExpression<VectorExpressionType::Add, Vector<T, Length, maxLength>, Vector<T, Length, maxLength>>
            operator+(const Vector<T, Length, maxLength>& v1, const Vector<T, Length, maxLength>& v2);

    template<class T, size_t Length, size_t maxLength>
    inline VectorExpression<VectorExpressionType::Sub, Vector<T, Length, maxLength>, Vector<T, Length, maxLength>>
            operator-(const Vector<T, Length, maxLength>& v1, const Vector<T, Length, maxLength>& v2);

    template<class T, size_t Length, size_t maxLength>
    T operator*(const Vector<T, Length, maxLength>& v1, const Vector<T, Length, maxLength>& v2);

    template<class T, size_t Length, size_t maxLength>
    VectorExpression<VectorExpressionType::Add, Vector<T, Length, maxLength>, T> operator+(const Vector<T, Length, maxLength>& v, const T& s);

    template<class T, size_t Length, size_t maxLength>
    VectorExpression<VectorExpressionType::Sub, Vector<T, Length, maxLength>, T> operator-(const Vector<T, Length, maxLength>& v, const T& s);

    template<class T, size_t Length, size_t maxLength>
    VectorExpression<VectorExpressionType::Mul, Vector<T, Length, maxLength>, T> operator*(const Vector<T, Length, maxLength>& v, const T& s);

    template<class T, size_t Length, size_t maxLength>
    VectorExpression<VectorExpressionType::Div, Vector<T, Length, maxLength>, T> operator/(const Vector<T, Length, maxLength>& v, const T& s);
    
    template<class T, size_t Length, size_t maxLength, VectorExpressionType type, class T1, class T2>
    void operator+=(Vector<T, Length, maxLength>& v1, const VectorExpression<type, T1, T2>& exp);

    template<class T, size_t Length, size_t maxLength, VectorExpressionType type, class T1, class T2>
    void operator-=(Vector<T, Length, maxLength>& v1, const VectorExpression<type, T1, T2>& exp);
    /* Inline Implements */
    template<class T, size_t Length, size_t maxLength>
    inline void operator+=(Vector<T, Length, maxLength>& v1, const Vector<T, Length, maxLength>& v2) { v1 = v1 + v2; }

    template<class T, size_t Length, size_t maxLength>
    inline void operator-=(Vector<T, Length, maxLength>& v1, const Vector<T, Length, maxLength>& v2) { v1 = v1 - v2; }

    template<class T, size_t Length, size_t maxLength>
    inline void operator+=(Vector<T, Length, maxLength>& v, const T& n) { v = v + n; }

    template<class T, size_t Length, size_t maxLength>
    inline void operator-=(Vector<T, Length, maxLength>& v, const T& n) { v = v - n; }

    template<class T, size_t Length, size_t maxLength>
    inline void operator*=(Vector<T, Length, maxLength>& v, const T& n) { v = v * n; }

    template<class T, size_t Length, size_t maxLength>
    inline void operator/=(Vector<T, Length, maxLength>& v, const T& n) { v = v / n; }
}

#include "VectorImpl.h"
