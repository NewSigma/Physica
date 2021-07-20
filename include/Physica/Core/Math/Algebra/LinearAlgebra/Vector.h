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
#include "VectorBlock.h"
#include "VectorExpression.h"
#include "Matrix/DenseMatrixImpl/DenseMatrixType.h"

namespace Physica::Core {
    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrix;
    /**
     * T must be either Scalar or ComplexScalar.
     * 
     * Default template arguments are defined in \file VectorExpression.h
     */
    template<class T, size_t Length, size_t MaxLength>
    class Vector : public Utils::Array<T, Length, MaxLength> {
    public:
        using ScalarType = T;
        using VectorType = Vector<T, Length, MaxLength>; //Redeclare self for the implementation of VectorExpression
    private:
        static_assert(Length == Dynamic || Length == MaxLength, "MaxLength of fixed vector must equals to its length.");
        using Base = Utils::Array<T, Length, MaxLength>;
        using ColMatrix = DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, Length, 1, MaxLength, 1>;
        using RowMatrix = DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, 1, Length, 1, MaxLength>;
    public:
        using Base::Base;
        Vector() = default;
        template<Utils::ExpressionType type, class T1, class T2>
        Vector(const VectorExpression<type, T1, T2>& expression); //NOLINT Implicit conversions is permitted.
        Vector(const Vector&) = default;
        Vector(Vector&&) noexcept = default;
        ~Vector() = default;
        /* Operators */
        Vector& operator=(const Vector&) = default;
        Vector& operator=(Vector&&) noexcept = default;
        template<Utils::ExpressionType type, class T1, class T2>
        Vector& operator=(const VectorExpression<type, T1, T2>& exp);
        /* Operations */
        Vector& toOpposite();
        void toUnit();
        ColMatrix copyToColMatrix() const;
        ColMatrix moveToColMatrix();
        RowMatrix copyToRowMatrix() const;
        RowMatrix moveToRowMatrix();
        /* Getters */
        [[nodiscard]] bool isZero() const;
        [[nodiscard]] T max() const;
        [[nodiscard]] T min() const;
        [[nodiscard]] T norm() const;
        [[nodiscard]] T squaredNorm() const;
        /* Helpers */
        static Vector<T> zeroVector(size_t len);
        static Vector<T> randomVector(size_t len);
        static Vector simplyMultiply(const Vector& v1, const Vector& v2);

        template<class Derived>
        friend class Internal::VectorExpressionHelper;
    };
    /* Operators */
    template<class T, size_t Length, size_t MaxLength>
    std::ostream& operator<<(std::ostream& os, const Vector<T, Length, MaxLength>& v);

    template<class T, size_t Length, size_t MaxLength>
    inline VectorExpression<Utils::ExpressionType::Minus, Vector<T, Length, MaxLength>> operator-(const Vector<T, Length, MaxLength>& v);

    template<class T, size_t Length, size_t MaxLength>
    inline VectorExpression<Utils::ExpressionType::Add, Vector<T, Length, MaxLength>, Vector<T, Length, MaxLength>>
            operator+(const Vector<T, Length, MaxLength>& v1, const Vector<T, Length, MaxLength>& v2);

    template<class T, size_t Length, size_t MaxLength>
    inline VectorExpression<Utils::ExpressionType::Sub, Vector<T, Length, MaxLength>, Vector<T, Length, MaxLength>>
            operator-(const Vector<T, Length, MaxLength>& v1, const Vector<T, Length, MaxLength>& v2);

    template<class T, size_t Length, size_t MaxLength>
    T operator*(const Vector<T, Length, MaxLength>& v1, const Vector<T, Length, MaxLength>& v2);

    template<class T, size_t Length, size_t MaxLength>
    VectorExpression<Utils::ExpressionType::Add, Vector<T, Length, MaxLength>, T> operator+(const Vector<T, Length, MaxLength>& v, const T& s);

    template<class T, size_t Length, size_t MaxLength>
    inline VectorExpression<Utils::ExpressionType::Add, Vector<T, Length, MaxLength>, T> operator+(const T& s, const Vector<T, Length, MaxLength>& v) { return v + s; }

    template<class T, size_t Length, size_t MaxLength>
    VectorExpression<Utils::ExpressionType::Sub, Vector<T, Length, MaxLength>, T> operator-(const Vector<T, Length, MaxLength>& v, const T& s);

    template<class T, size_t Length, size_t MaxLength>
    VectorExpression<Utils::ExpressionType::Mul, Vector<T, Length, MaxLength>, T> operator*(const Vector<T, Length, MaxLength>& v, const T& s);

    template<class T, size_t Length, size_t MaxLength>
    VectorExpression<Utils::ExpressionType::Mul, Vector<T, Length, MaxLength>, T> operator*(const T& s, const Vector<T, Length, MaxLength>& v) { return v * s; }

    template<class T, size_t Length, size_t MaxLength>
    VectorExpression<Utils::ExpressionType::Div, Vector<T, Length, MaxLength>, T> operator/(const Vector<T, Length, MaxLength>& v, const T& s);
    
    template<class T, size_t Length, size_t MaxLength, Utils::ExpressionType type, class T1, class T2>
    void operator+=(Vector<T, Length, MaxLength>& v1, const VectorExpression<type, T1, T2>& exp);

    template<class T, size_t Length, size_t MaxLength, Utils::ExpressionType type, class T1, class T2>
    void operator-=(Vector<T, Length, MaxLength>& v1, const VectorExpression<type, T1, T2>& exp);
    /* Inline Implements */
    template<class T, size_t Length, size_t MaxLength>
    inline void operator+=(Vector<T, Length, MaxLength>& v1, const Vector<T, Length, MaxLength>& v2) { v1 = v1 + v2; }

    template<class T, size_t Length, size_t MaxLength>
    inline void operator-=(Vector<T, Length, MaxLength>& v1, const Vector<T, Length, MaxLength>& v2) { v1 = v1 - v2; }

    template<class T, size_t Length, size_t MaxLength>
    inline void operator+=(Vector<T, Length, MaxLength>& v, const T& n) { v = v + n; }

    template<class T, size_t Length, size_t MaxLength>
    inline void operator-=(Vector<T, Length, MaxLength>& v, const T& n) { v = v - n; }

    template<class T, size_t Length, size_t MaxLength>
    inline void operator*=(Vector<T, Length, MaxLength>& v, const T& n) { v = v * n; }

    template<class T, size_t Length, size_t MaxLength>
    inline void operator/=(Vector<T, Length, MaxLength>& v, const T& n) { v = v / n; }
}

#include "VectorImpl.h"
