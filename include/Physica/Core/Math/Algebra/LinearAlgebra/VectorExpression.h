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
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Utils/Template/CRTPBase.h"
/**
 * This file contains implementation of expression templates of \class Vector.
 * 
 * Optimize: Add compile time length check, length of two vectors must be equal.
 */
namespace Physica::Core {
    //Forward declaration
    using Utils::Dynamic;
    template<class T = MultiScalar, size_t Length = Dynamic, size_t MaxLength = Length>
    class Vector;

    enum class VectorExpressionType {
        Minus,
        Add,
        Sub,
        Mul,
        Div
    };
    /**
     * \class VectorExpression represents \param T1 \param type \param T2. e.g. vector + scalar, expression * expression
     */
    template<VectorExpressionType type, class T1, class T2 = T1>
    class VectorExpression {
        //Here the condition must always be false.
        static_assert(type > VectorExpressionType::Div, "Not implemented.");
    };

    namespace Internal {
        template<class T>
        class Traits;
        /**
         * VectorType: Type of the vector calculated from the expression.
         */
        template<VectorExpressionType type, class T, size_t Length, size_t MaxLength>
        class Traits<VectorExpression<type, Vector<T, Length, MaxLength>, Vector<T, Length, MaxLength>>> {
        public:
            using VectorType = Vector<T, Length, MaxLength>;
        };

        template<VectorExpressionType type, class T, size_t Length, size_t MaxLength, ScalarType scalarType, bool errorTrack>
        class Traits<VectorExpression<type, Vector<T, Length, MaxLength>, Scalar<scalarType, errorTrack>>> {
        public:
            using VectorType = Vector<T, Length, MaxLength>;
        };

        template<VectorExpressionType type, VectorExpressionType type1, class T11, class T12
                                            , VectorExpressionType type2, class T21, class T22>
        class Traits<VectorExpression<type, VectorExpression<type1, T11, T12>, VectorExpression<type2, T21, T22>>> {
            static_assert(std::is_same<typename VectorExpression<type1, T11, T12>::VectorType
                                       , typename VectorExpression<type2, T21, T22>::VectorType>::value
                          , "Types of the two operands must be same.");
        public:
            using VectorType = typename VectorExpression<type1, T11, T12>::VectorType;
        };

        template<VectorExpressionType type, VectorExpressionType type1, class T1, class T2, class T, size_t Length, size_t MaxLength>
        class Traits<VectorExpression<type, VectorExpression<type1, T1, T2>, Vector<T, Length, MaxLength>>> {
            static_assert(std::is_same<typename VectorExpression<type1, T1, T2>::VectorType, Vector<T, Length, MaxLength>>::value
                          , "Types of the two operands must be same.");
        public:
            using VectorType = Vector<T, Length, MaxLength>;
        };

        template<VectorExpressionType type, VectorExpressionType type1, class T1, class T2, class T, size_t Length, size_t MaxLength>
        class Traits<VectorExpression<type, Vector<T, Length, MaxLength>, VectorExpression<type1, T1, T2>>> {
            static_assert(std::is_same<typename VectorExpression<type1, T1, T2>::VectorType, Vector<T, Length, MaxLength>>::value
                          , "Types of the two operands must be same.");
        public:
            using VectorType = Vector<T, Length, MaxLength>;
        };

        template<VectorExpressionType type, VectorExpressionType type1, class T1, class T2, ScalarType scalarType, bool errorTrack>
        class Traits<VectorExpression<type, VectorExpression<type1, T1, T2>, Scalar<scalarType, errorTrack>>> {
        public:
            using VectorType = typename VectorExpression<type1, T1, T2>::VectorType;
        };
        /**
         * This class implements calc() for all \class VectorExpression
         */
        template<class Derived>
        class VectorExpressionHelper : public Utils::CRTPBase<Derived> {
        public:
            using VectorType = typename Traits<Derived>::VectorType;
            using ScalarType = typename VectorType::ScalarType;
        private:
            using Base = Utils::CRTPBase<Derived>;
        public:
            /* Operators */
            operator VectorType() { return calc(); }
            /* Getters */
            [[nodiscard]] VectorType calc() const {
                VectorType result{};
                const Derived& derived = Base::getDerived();
                const size_t length = derived.getLength();
                result.reserve(length);
                for (size_t i = 0; i < length; ++i)
                    result.allocate(derived[i], i);
                result.setLength(length);
                return result;
            }
        };
    }
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class T>
    class VectorExpression<VectorExpressionType::Minus, T>
            : public Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Minus, T>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Minus, T>>;
        const T& exp;
    public:
        explicit VectorExpression(const T& exp_) : exp(exp_) {}

        typename Base::ScalarType operator[](size_t s) const { return -exp[s]; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Add//////////////////////////////////////
    template<class T1, class T2>
    class VectorExpression<VectorExpressionType::Add, T1, T2>
            : public Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Add, T1, T2>> {
        static_assert(std::is_same<typename T1::VectorType, typename T2::VectorType>::value
                      , "Types of two operands of add must be same.");
        using Base = Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Add, T1, T2>>;
        const T1& exp1;
        const T2& exp2;
    public:
        VectorExpression(const T1& exp1_, const T2& exp2_) : exp1(exp1_), exp2(exp2_) {
            assert(exp1.getLength() == exp2.getLength());
        }

        typename Base::ScalarType operator[](size_t s) const { return exp1[s] + exp2[s]; }
        [[nodiscard]] size_t getLength() const { return exp1.getLength(); }
    };

    template<class T, ScalarType type, bool errorTrack>
    class VectorExpression<VectorExpressionType::Add, T, Scalar<type, errorTrack>>
            : public Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Add, T, Scalar<type, errorTrack>>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Add, T, Scalar<type, errorTrack>>>;
        const T& exp;
        const Scalar<type, errorTrack>& scalar;
    public:
        VectorExpression(const T& exp_, const Scalar<type, errorTrack>& scalar_) : exp(exp_), scalar(scalar_) {}

        typename Base::ScalarType operator[](size_t s) const { return exp[s] + scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class T1, class T2>
    class VectorExpression<VectorExpressionType::Sub, T1, T2>
            : public Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Sub, T1, T2>> {
        static_assert(std::is_same<typename T1::VectorType, typename T2::VectorType>::value
                      , "Types of two operands of sub must be same.");
        using Base = Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Sub, T1, T2>>;
        T1 exp1;
        T2 exp2;
    public:
        VectorExpression(const T1& exp1_, const T2& exp2_) : exp1(exp1_), exp2(exp2_) {
            assert(exp1.getLength() == exp2.getLength());
        }

        typename Base::ScalarType operator[](size_t s) const { return exp1[s] - exp2[s]; }
        [[nodiscard]] size_t getLength() const { return exp1.getLength(); }
    };

    template<class T, ScalarType scalarType, bool errorTrack>
    class VectorExpression<VectorExpressionType::Sub, T, Scalar<scalarType, errorTrack>>
            : public Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Sub, T, Scalar<scalarType, errorTrack>>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Sub, T, Scalar<scalarType, errorTrack>>>;
        T exp;
        const Scalar<scalarType, errorTrack>& scalar;
    public:
        VectorExpression(const T& exp_, const Scalar<scalarType, errorTrack>& scalar_) : exp(exp_), scalar(scalar_) {}

        typename Base::ScalarType operator[](size_t s) const { return exp[s] - scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class T, ScalarType scalarType, bool errorTrack>
    class VectorExpression<VectorExpressionType::Mul, T, Scalar<scalarType, errorTrack>>
            : public Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Mul, T, Scalar<scalarType, errorTrack>>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Mul, T, Scalar<scalarType, errorTrack>>>;
        T exp;
        const Scalar<scalarType, errorTrack>& scalar;
    public:
        VectorExpression(const T& exp_, const Scalar<scalarType, errorTrack>& scalar_) : exp(exp_), scalar(scalar_) {}

        typename Base::ScalarType operator[](size_t s) const { return exp[s] * scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Div//////////////////////////////////////
    template<class T, ScalarType scalarType, bool errorTrack>
    class VectorExpression<VectorExpressionType::Div, T, Scalar<scalarType, errorTrack>>
            : public Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Div, T, Scalar<scalarType, errorTrack>>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<VectorExpressionType::Div, T, Scalar<scalarType, errorTrack>>>;
        T exp;
        const Scalar<scalarType, errorTrack>& scalar;
    public:
        VectorExpression(const T& exp_, const Scalar<scalarType, errorTrack>& scalar_) : exp(exp_), scalar(scalar_) {}

        typename Base::ScalarType operator[](size_t s) const { return exp[s] / scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    template<VectorExpressionType type, class T1, class T2>
    inline VectorExpression<VectorExpressionType::Minus, VectorExpression<type, T1, T2>>
    operator-(const VectorExpression<type, T1, T2>& exp) {
        return VectorExpression<VectorExpressionType::Minus
                                , VectorExpression<type, T1, T2>>(exp);
    }
    //////////////////////////////////////Add//////////////////////////////////////
    template<VectorExpressionType type, class T1, class T2, ScalarType scalarType, bool errorTrack>
    inline VectorExpression<VectorExpressionType::Add, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>>
    operator+(const VectorExpression<type, T1, T2>& exp, const Scalar<scalarType, errorTrack>& s) {
        return VectorExpression<VectorExpressionType::Add
                                , VectorExpression<type, T1, T2>
                                , Scalar<scalarType, errorTrack>>(exp, s);
    }

    template<VectorExpressionType type, class T1, class T2, class T, size_t Length, size_t MaxLength>
    inline VectorExpression<VectorExpressionType::Add, VectorExpression<type, T1, T2>, Vector<T, Length, MaxLength>>
    operator+(const VectorExpression<type, T1, T2>& exp, const Vector<T, Length, MaxLength>& v) {
        return VectorExpression<VectorExpressionType::Add
                                , VectorExpression<type, T1, T2>
                                , Vector<T, Length, MaxLength>>(exp, v);
    }

    template<VectorExpressionType type, class T1, class T2, class T, size_t Length, size_t MaxLength>
    inline VectorExpression<VectorExpressionType::Add, VectorExpression<type, T1, T2>, Vector<T, Length, MaxLength>>
    operator+(const Vector<T, Length, MaxLength>& v, const VectorExpression<type, T1, T2>& exp) {
        return exp + v;
    }

    template<VectorExpressionType type1, class T11, class T12, VectorExpressionType type2, class T21, class T22>
    inline VectorExpression<VectorExpressionType::Add, VectorExpression<type1, T11, T12>, VectorExpression<type2, T21, T22>>
    operator+(const VectorExpression<type1, T11, T12>& exp1, const VectorExpression<type2, T21, T22>& exp2) {
        return VectorExpression<VectorExpressionType::Add
                                , VectorExpression<type1, T11, T12>
                                , VectorExpression<type2, T21, T22>>(exp1, exp2);
    }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<VectorExpressionType type, class T1, class T2, ScalarType scalarType, bool errorTrack>
    inline VectorExpression<VectorExpressionType::Sub, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>>
    operator-(const VectorExpression<type, T1, T2>& exp, const Scalar<scalarType, errorTrack>& s) {
        return VectorExpression<VectorExpressionType::Sub
                                , VectorExpression<type, T1, T2>
                                , Scalar<scalarType, errorTrack>>(exp, s);
    }

    template<VectorExpressionType type, class T1, class T2, class T, size_t Length, size_t MaxLength>
    inline VectorExpression<VectorExpressionType::Sub, VectorExpression<type, T1, T2>, Vector<T, Length, MaxLength>>
    operator-(const VectorExpression<type, T1, T2>& exp, const Vector<T, Length, MaxLength>& v) {
        return VectorExpression<VectorExpressionType::Sub
                                , VectorExpression<type, T1, T2>
                                , Vector<T, Length, MaxLength>>(exp, v);
    }

    template<VectorExpressionType type, class T1, class T2, class T, size_t Length, size_t MaxLength>
    inline VectorExpression<VectorExpressionType::Sub, Vector<T, Length, MaxLength>, VectorExpression<type, T1, T2>>
    operator-(const Vector<T, Length, MaxLength>& v, const VectorExpression<type, T1, T2>& exp) {
        return VectorExpression<VectorExpressionType::Sub
                                , Vector<T, Length, MaxLength>
                                , VectorExpression<type, T1, T2>>(v, exp);
    }

    template<VectorExpressionType type1, class T11, class T12, VectorExpressionType type2, class T21, class T22>
    inline VectorExpression<VectorExpressionType::Sub, VectorExpression<type1, T11, T12>, VectorExpression<type2, T21, T22>>
    operator-(const VectorExpression<type1, T11, T12>& exp1, const VectorExpression<type2, T21, T22>& exp2) {
        return VectorExpression<VectorExpressionType::Sub
                                , VectorExpression<type1, T11, T12>
                                , VectorExpression<type2, T21, T22>>(exp1, exp2);
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<VectorExpressionType type, class T1, class T2, ScalarType scalarType, bool errorTrack>
    inline VectorExpression<VectorExpressionType::Mul, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>>
    operator*(const VectorExpression<type, T1, T2>& exp, const Scalar<scalarType, errorTrack>& s) {
        return VectorExpression<VectorExpressionType::Mul
                                , VectorExpression<type, T1, T2>
                                , Scalar<scalarType, errorTrack>>(exp, s);
    }
    //////////////////////////////////////Div//////////////////////////////////////
    template<VectorExpressionType type, class T1, class T2, ScalarType scalarType, bool errorTrack>
    inline VectorExpression<VectorExpressionType::Div, VectorExpression<type, T1, T2>, Scalar<scalarType, errorTrack>>
    operator/(const VectorExpression<type, T1, T2>& exp, const Scalar<scalarType, errorTrack>& s) {
        return VectorExpression<VectorExpressionType::Div
                                , VectorExpression<type, T1, T2>
                                , Scalar<scalarType, errorTrack>>(exp, s);
    }
}
