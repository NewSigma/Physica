/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Utils/Template/ExpressionTemplateHelper.h"

namespace Physica::Core {
    /**
     * \class DenseMatrixExpression represents \param T1 \param type \param T2. e.g. matrix + scalar, expression * expression
     */
    template<Utils::ExpressionType type, class T1, class T2 = T1>
    class DenseMatrixExpression;
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class T>
    class DenseMatrixExpression<Utils::ExpressionType::Minus, T> {
    public:
        using ScalarType = typename T::ScalarType;
    private:
        const T& exp;
    public:
        DenseMatrixExpression(const T& exp_) : exp(exp_) {}

        [[nodiscard]] ScalarType operator()(size_t row, size_t col) const { return -exp(row, col); }
        [[nodiscard]] size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Add//////////////////////////////////////
    template<class T1, class T2>
    class DenseMatrixExpression<Utils::ExpressionType::Add, T1, T2> {
    public:
        using ScalarType = typename Internal::BinaryScalarOpReturnType<typename T1::ScalarType, typename T2::ScalarType>::Type;
    private:
        const T1& exp1;
        const T2& exp2;
    public:
        DenseMatrixExpression(const T1& exp1_, const T2& exp2_) : exp1(exp1_), exp2(exp2_) {}

        [[nodiscard]] ScalarType operator()(size_t row, size_t col) const { return ScalarType(exp1(row, col)) + ScalarType(exp2(row, col)); }
        [[nodiscard]] size_t getRow() const { return exp1.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp1.getColumn(); }
    };

    template<class T, class AnyScalar>
    class DenseMatrixExpression<Utils::ExpressionType::Add, T, ScalarBase<AnyScalar>> {
    public:
        using ScalarType = typename Internal::BinaryScalarOpReturnType<typename T::ScalarType, AnyScalar>::Type;
    private:
        const T& exp;
        const AnyScalar& scalar;
    public:
        DenseMatrixExpression(const T& exp_, const ScalarBase<AnyScalar>& base) : exp(exp_), scalar(base.getDerived()) {}

        [[nodiscard]] ScalarType operator()(size_t row, size_t col) const { return ScalarType(exp(row, col)) + ScalarType(scalar); }
        [[nodiscard]] size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class T1, class T2>
    class DenseMatrixExpression<Utils::ExpressionType::Sub, T1, T2> {
    public:
        using ScalarType = typename Internal::BinaryScalarOpReturnType<typename T1::ScalarType, typename T2::ScalarType>::Type;
    private:
        const T1& exp1;
        const T2& exp2;
    public:
        DenseMatrixExpression(const T1& exp1_, const T2& exp2_) : exp1(exp1_), exp2(exp2_) {}

        [[nodiscard]] ScalarType operator()(size_t row, size_t col) const { return ScalarType(exp1(row, col)) - ScalarType(exp2(row, col)); }
        [[nodiscard]] size_t getRow() const { return exp1.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp1.getColumn(); }
    };

    template<class T, class AnyScalar>
    class DenseMatrixExpression<Utils::ExpressionType::Sub, T, ScalarBase<AnyScalar>> {
    public:
        using ScalarType = typename Internal::BinaryScalarOpReturnType<typename T::ScalarType, AnyScalar>::Type;
    private:
        const T& exp;
        const AnyScalar& scalar;
    public:
        DenseMatrixExpression(const T& exp_, const ScalarBase<AnyScalar>& base) : exp(exp_), scalar(base.getDerived()) {}

        [[nodiscard]] ScalarType operator()(size_t row, size_t col) const { return ScalarType(exp(row, col)) - ScalarType(scalar); }
        [[nodiscard]] size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class T, class AnyScalar>
    class DenseMatrixExpression<Utils::ExpressionType::Mul, T, ScalarBase<AnyScalar>> {
    public:
        using ScalarType = typename Internal::BinaryScalarOpReturnType<typename T::ScalarType, AnyScalar>::Type;
    private:
        const T& exp;
        const AnyScalar& scalar;
    public:
        DenseMatrixExpression(const T& exp_, const ScalarBase<AnyScalar>& base) : exp(exp_), scalar(base.getDerived()) {}

        [[nodiscard]] ScalarType operator()(size_t row, size_t col) const { return ScalarType(exp(row, col)) * ScalarType(scalar); }
        [[nodiscard]] size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Div//////////////////////////////////////
    template<class T, class AnyScalar>
    class DenseMatrixExpression<Utils::ExpressionType::Div, T, ScalarBase<AnyScalar>> {
    public:
        using ScalarType = typename Internal::BinaryScalarOpReturnType<typename T::ScalarType, AnyScalar>::Type;
    private:
        const T& exp;
        const AnyScalar& scalar;
    public:
        DenseMatrixExpression(const T& exp_, const ScalarBase<AnyScalar>& base) : exp(exp_), scalar(base.getDerived()) {}

        [[nodiscard]] ScalarType operator()(size_t row, size_t col) const { return exp(row, col) / ScalarType(scalar); }
        [[nodiscard]] size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] size_t getColumn() const { return exp.getColumn(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class Derived>
    inline DenseMatrixExpression<Utils::ExpressionType::Minus, DenseMatrixBase<Derived>>
    operator-(const DenseMatrixBase<Derived>& mat) {
        return DenseMatrixExpression<Utils::ExpressionType::Minus, DenseMatrixBase<Derived>>(mat);
    }

    template<Utils::ExpressionType type, class T1, class T2>
    inline DenseMatrixExpression<Utils::ExpressionType::Minus, DenseMatrixExpression<type, T1, T2>>
    operator-(const DenseMatrixExpression<type, T1, T2>& exp) {
        return DenseMatrixExpression<Utils::ExpressionType::Minus, DenseMatrixExpression<type, T1, T2>>(exp);
    }
    //////////////////////////////////////Add//////////////////////////////////////
    template<class Derived, class Exp>
    inline DenseMatrixExpression<Utils::ExpressionType::Add, DenseMatrixBase<Derived>, Exp>
    operator+(const DenseMatrixBase<Derived>& mat, const Exp& exp) {
        return DenseMatrixExpression<Utils::ExpressionType::Add, DenseMatrixBase<Derived>, Exp>(mat, exp);
    }

    template<Utils::ExpressionType type, class T1, class T2, class Exp>
    inline DenseMatrixExpression<Utils::ExpressionType::Add, DenseMatrixExpression<type, T1, T2>, Exp>
    operator+(const DenseMatrixExpression<type, T1, T2>& exp1, const Exp& exp2) {
        return DenseMatrixExpression<Utils::ExpressionType::Add, DenseMatrixExpression<type, T1, T2>, Exp>(exp1, exp2);
    }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class Derived, class Exp>
    inline DenseMatrixExpression<Utils::ExpressionType::Sub, DenseMatrixBase<Derived>, Exp>
    operator-(const DenseMatrixBase<Derived>& mat, const Exp& exp) {
        return DenseMatrixExpression<Utils::ExpressionType::Sub, DenseMatrixBase<Derived>, Exp>(mat, exp);
    }

    template<Utils::ExpressionType type, class T1, class T2, class Exp>
    inline DenseMatrixExpression<Utils::ExpressionType::Sub, DenseMatrixExpression<type, T1, T2>, Exp>
    operator-(const DenseMatrixExpression<type, T1, T2>& exp1, const Exp& exp2) {
        return DenseMatrixExpression<Utils::ExpressionType::Sub, DenseMatrixExpression<type, T1, T2>, Exp>(exp1, exp2);
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<Utils::ExpressionType type, class T1, class T2, class Exp>
    inline DenseMatrixExpression<Utils::ExpressionType::Mul, DenseMatrixExpression<type, T1, T2>, Exp>
    operator*(const DenseMatrixExpression<type, T1, T2>& exp1, const Exp& exp2) {
        return DenseMatrixExpression<Utils::ExpressionType::Mul, DenseMatrixExpression<type, T1, T2>, Exp>(exp1, exp2);
    }
    //////////////////////////////////////Div//////////////////////////////////////
    template<Utils::ExpressionType type, class T1, class T2, class Exp>
    inline DenseMatrixExpression<Utils::ExpressionType::Div, DenseMatrixExpression<type, T1, T2>, Exp>
    operator/(const DenseMatrixExpression<type, T1, T2>& exp1, const Exp& exp2) {
        return DenseMatrixExpression<Utils::ExpressionType::Div, DenseMatrixExpression<type, T1, T2>, Exp>(exp1, exp2);
    }
}
