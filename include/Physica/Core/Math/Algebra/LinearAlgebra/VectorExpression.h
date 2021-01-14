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
#ifndef PHYSICA_VECTOREXPRESSION_H
#define PHYSICA_VECTOREXPRESSION_H

#include "Physica/Core/MultiPrecision/Scalar.h"
/**
 * This file contains implementation of expression templates of @class Vector.
 */
namespace Physica::Core {
    template<class T = MultiScalar, size_t maxLength = Utils::Dynamic>
    class Vector;

    enum class VectorExpressionType {
        Minus,
        Add,
        Sub,
        Mul,
        Div
    };
    /*!
     * @class VectorExpression represents @param T1 @param type @param T2. e.g. vector + scalar, expression * expression
     */
    template<VectorExpressionType type, class T1, class T2 = T1>
    class VectorExpression {
        //Here the condition must always be false.
        static_assert(type > VectorExpressionType::Div, "Not implemented.");
    };
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class T, size_t maxLength>
    class VectorExpression<VectorExpressionType::Minus, Vector<T, maxLength>> {
        const Vector<T, maxLength>& vector;
    public:
        explicit VectorExpression(const Vector<T, maxLength>& v) : vector(v) {}

        const T& operator[](size_t s) const { return vector[s]; }
        [[nodiscard]] size_t getLength() const { return vector.getLength(); }
    };

    template<VectorExpressionType type, class T1, class T2>
    class VectorExpression<VectorExpressionType::Minus, VectorExpression<type, T1, T2>> {
        VectorExpression<type, T1, T2> exp;
    public:
        explicit VectorExpression(const VectorExpression<type, T1, T2>& exp) : exp(exp) {}

        decltype(auto) operator[](size_t s) { return -exp[s]; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Add//////////////////////////////////////
    template<class T, size_t maxLength>
    class VectorExpression<VectorExpressionType::Add, Vector<T, maxLength>, Vector<T, maxLength>> {
        const Vector<T, maxLength>& v1;
        const Vector<T, maxLength>& v2;
    public:
        explicit VectorExpression(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2) : v1(v1), v2(v2) {}

        decltype(auto) operator[](size_t s) const { return v1[s] + v2[s]; }
        [[nodiscard]] size_t getLength() const { return v1.getLength(); }
    };

    template<class T, size_t maxLength, ScalarType type, bool errorTrack>
    class VectorExpression<VectorExpressionType::Add, Vector<T, maxLength>, Scalar < type, errorTrack>> {
        const Vector<T, maxLength>& v;
        const Scalar<type, errorTrack>& scalar;
    public:
        explicit VectorExpression(const Vector<T, maxLength>& v, const Scalar<type, errorTrack>& scalar) : v(v), scalar(scalar) {}

        decltype(auto) operator[](size_t s) const { return v[s] + scalar; }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<VectorExpressionType type, class T1, class T2, class T, size_t maxLength>
    class VectorExpression<VectorExpressionType::Add, VectorExpression<type, T1, T2>, Vector<T, maxLength>> {
        VectorExpression<type, T1, T2> exp;
        const Vector<T, maxLength>& v;
    public:
        explicit VectorExpression(const VectorExpression<type, T1, T2>& exp, const Vector<T, maxLength>& v) : exp(exp), v(v) {}

        decltype(auto) operator[](size_t s) const { return  exp[s] + v[s]; }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<ScalarType scalarType, bool errorTrack, VectorExpressionType type, class T1, class T2>
    class VectorExpression<VectorExpressionType::Add, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>> {
        VectorExpression<type, T1, T2> exp;
        const Scalar<scalarType, errorTrack>& scalar;
    public:
        explicit VectorExpression(const VectorExpression<type, T1, T2>& exp, const Scalar<scalarType, errorTrack>& scalar) : exp(exp), scalar(scalar) {}

        decltype(auto) operator[](size_t s) const { return exp[s] + scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };

    template<VectorExpressionType type1, class T11, class T12, VectorExpressionType type2, class T21, class T22>
    class VectorExpression<VectorExpressionType::Add, VectorExpression<type1, T11, T12>, VectorExpression<type2, T21, T22>> {
        VectorExpression<type1, T11, T12> exp1;
        VectorExpression<type2, T21, T22> exp2;
    public:
        explicit VectorExpression(const VectorExpression<type1, T11, T12>& exp1, const VectorExpression<type2, T21, T22>& exp2) : exp1(exp1), exp2(exp2) {}

        decltype(auto) operator[](size_t s) const { return exp1[s] + exp2[s]; }
        [[nodiscard]] size_t getLength() const { return exp1.getLength(); }
    };
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class T, size_t maxLength>
    class VectorExpression<VectorExpressionType::Sub, Vector<T, maxLength>, Vector<T, maxLength>> {
        const Vector<T, maxLength>& v1;
        const Vector<T, maxLength>& v2;
    public:
        explicit VectorExpression(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2) : v1(v1), v2(v2) {}

        decltype(auto) operator[](size_t s) const { return v1[s] - v2[s]; }
        [[nodiscard]] size_t getLength() const { return v1.getLength(); }
    };

    template<class T, size_t maxLength, ScalarType type, bool errorTrack>
    class VectorExpression<VectorExpressionType::Sub, Vector<T, maxLength>, Scalar<type, errorTrack>> {
        const Vector<T, maxLength>& v;
        const Scalar<type, errorTrack>& scalar;
    public:
        explicit VectorExpression(const Vector<T, maxLength>& v, const Scalar<type, errorTrack>& scalar) : v(v), scalar(scalar) {}

        decltype(auto) operator[](size_t s) const { return v[s] - scalar; }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class T, size_t maxLength, VectorExpressionType type, class T1, class T2>
    class VectorExpression<VectorExpressionType::Sub, VectorExpression<type, T1, T2>, Vector<T, maxLength>> {
        VectorExpression<type, T1, T2> exp;
        const Vector<T, maxLength>& v;
    public:
        explicit VectorExpression(const VectorExpression<type, T1, T2>& exp, const Vector<T, maxLength>& v) : exp(exp), v(v) {}

        decltype(auto) operator[](size_t s) const { return exp[s] - v[s]; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };

    template<VectorExpressionType type, class T1, class T2, ScalarType scalarType, bool errorTrack>
    class VectorExpression<VectorExpressionType::Sub, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>> {
        VectorExpression<type, T1, T2> exp;
        const Scalar<scalarType, errorTrack>& scalar;
    public:
        explicit VectorExpression(const VectorExpression<type, T1, T2>& exp, const Scalar<scalarType, errorTrack>& scalar) : exp(exp), scalar(scalar) {}

        decltype(auto) operator[](size_t s) const { return exp[s] - scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };

    template<VectorExpressionType type1, class T11, class T12, VectorExpressionType type2, class T21, class T22>
    class VectorExpression<VectorExpressionType::Sub, VectorExpression<type1, T11, T12>, VectorExpression<type2, T21, T22>> {
        VectorExpression<type1, T11, T12> exp1;
        VectorExpression<type2, T21, T22> exp2;
    public:
        explicit VectorExpression(const VectorExpression<type1, T11, T12>& exp1, const VectorExpression<type2, T21, T22>& exp2) : exp1(exp1), exp2(exp2) {}

        decltype(auto) operator[](size_t s) const { return exp1[s] - exp2[s]; }
        [[nodiscard]] size_t getLength() const { return exp1.getLength(); }
    };
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class T, size_t maxLength, ScalarType type, bool errorTrack>
    class VectorExpression<VectorExpressionType::Mul, Vector<T, maxLength>, Scalar<type, errorTrack>> {
        const Vector<T, maxLength>& v;
        const Scalar<type, errorTrack>& scalar;
    public:
        explicit VectorExpression(const Vector<T, maxLength>& v, const Scalar<type, errorTrack>& scalar) : v(v), scalar(scalar) {}

        decltype(auto) operator[](size_t s) const { return v[s] * scalar; }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<ScalarType scalarType, bool errorTrack, VectorExpressionType type, class T1, class T2>
    class VectorExpression<VectorExpressionType::Mul, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>> {
        VectorExpression<type, T1, T2> exp;
        const Scalar<scalarType, errorTrack>& scalar;
    public:
        explicit VectorExpression(const VectorExpression<type, T1, T2>& exp, const Scalar<scalarType, errorTrack>& scalar) : exp(exp), scalar(scalar) {}

        decltype(auto) operator[](size_t s) const { return exp[s] * scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Div//////////////////////////////////////
    template<class T, size_t maxLength, ScalarType type, bool errorTrack>
    class VectorExpression<VectorExpressionType::Div, Vector<T, maxLength>, Scalar<type, errorTrack>> {
        const Vector<T, maxLength>& v;
        const Scalar<type, errorTrack>& scalar;
    public:
        explicit VectorExpression(const Vector<T, maxLength>& v, const Scalar<type, errorTrack>& scalar) : v(v), scalar(scalar) {}

        decltype(auto) operator[](size_t s) const { return v[s] / scalar; }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<ScalarType scalarType, bool errorTrack, VectorExpressionType type, class T1, class T2>
    class VectorExpression<VectorExpressionType::Div, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>> {
        VectorExpression<type, T1, T2> exp;
        const Scalar<scalarType, errorTrack>& scalar;
    public:
        explicit VectorExpression(const VectorExpression<type, T1, T2>& exp, const Scalar<scalarType, errorTrack>& scalar) : exp(exp), scalar(scalar) {}

        decltype(auto) operator[](size_t s) const { return exp[s] / scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    template<VectorExpressionType type, class T1, class T2>
    VectorExpression<VectorExpressionType::Minus, VectorExpression<type, T1, T2>>
    operator-(const VectorExpression<type, T1, T2>& exp) {
        return VectorExpression<VectorExpressionType::Minus, VectorExpression<type, T1, T2>>(exp);
    }
    //////////////////////////////////////Add//////////////////////////////////////
    template<VectorExpressionType type, class T1, class T2, ScalarType scalarType, bool errorTrack>
    VectorExpression<VectorExpressionType::Add, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>>
    operator+(const VectorExpression<type, T1, T2>& exp, const Scalar<scalarType, errorTrack>& s) {
        return VectorExpression<VectorExpressionType::Add, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>>(exp, s);
    }

    template<VectorExpressionType type, class T1, class T2, class T, size_t maxLength>
    VectorExpression<VectorExpressionType::Add, VectorExpression<type, T1, T2>, Vector<T, maxLength>>
    operator+(const VectorExpression<type, T1, T2>& exp, const Vector<T, maxLength>& v) {
        return VectorExpression<VectorExpressionType::Add, VectorExpression<type, T1, T2>, Vector<T, maxLength>>(exp, v);
    }

    template<VectorExpressionType type1, class T11, class T12, VectorExpressionType type2, class T21, class T22>
    VectorExpression<VectorExpressionType::Add, VectorExpression<type1, T11, T12>, VectorExpression<type2, T21, T22>>
    operator+(const VectorExpression<type1, T11, T12>& exp1, const VectorExpression<type2, T21, T22>& exp2) {
        return VectorExpression<VectorExpressionType::Add, VectorExpression<type1, T11, T12>, VectorExpression<type2, T21, T22>>(exp1, exp2);
    }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<VectorExpressionType type, class T1, class T2, ScalarType scalarType, bool errorTrack>
    VectorExpression<VectorExpressionType::Sub, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>>
    operator-(const VectorExpression<type, T1, T2>& exp, const Scalar<scalarType, errorTrack>& s) {
        return VectorExpression<VectorExpressionType::Sub, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>>(exp, s);
    }

    template<VectorExpressionType type, class T1, class T2, class T, size_t maxLength>
    VectorExpression<VectorExpressionType::Sub, VectorExpression<type, T1, T2>, Vector<T, maxLength>>
    operator-(const VectorExpression<type, T1, T2>& exp, const Vector<T, maxLength>& v) {
        return VectorExpression<VectorExpressionType::Sub, VectorExpression<type, T1, T2>, Vector<T, maxLength>>(exp, v);
    }

    template<VectorExpressionType type1, class T11, class T12, VectorExpressionType type2, class T21, class T22>
    VectorExpression<VectorExpressionType::Sub, VectorExpression<type1, T11, T12>, VectorExpression<type2, T21, T22>>
    operator-(const VectorExpression<type1, T11, T12>& exp1, const VectorExpression<type2, T21, T22>& exp2) {
        return VectorExpression<VectorExpressionType::Sub, VectorExpression<type1, T11, T12>, VectorExpression<type2, T21, T22>>(exp1, exp2);
    }
}

#endif
