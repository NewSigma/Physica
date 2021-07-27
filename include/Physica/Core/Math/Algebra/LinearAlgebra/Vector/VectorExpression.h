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
#include "Physica/Utils/Template/ExpressionTemplateHelper.h"
#include "Physica/Utils/Template/CRTPBase.h"
/**
 * This file contains implementation of expression templates of \class Vector.
 * 
 * Optimize: Add compile time length check, length of two vectors must be equal.
 */
namespace Physica::Core {
    //Forward declaration
    using Utils::Dynamic;
    template<class T, size_t Length, size_t MaxLength>
    class Vector;
    /**
     * \class VectorExpression represents \param T1 \param type \param T2. e.g. vector + scalar, expression * expression
     */
    template<Utils::ExpressionType type, class T1, class T2 = T1>
    class VectorExpression;

    namespace Internal {
        template<class T>
        class Traits;
        /**
         * VectorType: Type of the vector calculated from the expression.
         */
        template<Utils::ExpressionType type, class Exp1, class Exp2>
        class Traits<VectorExpression<type, Exp1, Exp2>> {
            static_assert(std::is_same<typename Exp1::VectorType, typename Exp2::VectorType>::value
                          , "Types of the two operands must be same.");
        public:
            using VectorType = typename Exp1::VectorType;
        };

        template<Utils::ExpressionType type, class Exp, ScalarOption option, bool errorTrack>
        class Traits<VectorExpression<type, Exp, Scalar<option, errorTrack>>> {
        public:
            using VectorType = typename Exp::VectorType;
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
    class VectorExpression<Utils::ExpressionType::Minus, T>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Minus, T>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Minus, T>>;
        const T& exp;
    public:
        explicit VectorExpression(const T& exp_) : exp(exp_) {}

        typename Base::ScalarType operator[](size_t s) const { return -exp[s]; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Add//////////////////////////////////////
    template<class T1, class T2>
    class VectorExpression<Utils::ExpressionType::Add, T1, T2>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Add, T1, T2>> {
        static_assert(std::is_same<typename T1::VectorType, typename T2::VectorType>::value
                      , "Types of two operands of add must be same.");
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Add, T1, T2>>;
        const T1& exp1;
        const T2& exp2;
    public:
        VectorExpression(const T1& exp1_, const T2& exp2_) : exp1(exp1_), exp2(exp2_) {
            assert(exp1.getLength() == exp2.getLength());
        }

        typename Base::ScalarType operator[](size_t s) const { return exp1[s] + exp2[s]; }
        [[nodiscard]] size_t getLength() const { return exp1.getLength(); }
    };

    template<class T, ScalarOption option, bool errorTrack>
    class VectorExpression<Utils::ExpressionType::Add, T, Scalar<option, errorTrack>>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Add, T, Scalar<option, errorTrack>>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Add, T, Scalar<option, errorTrack>>>;
        const T& exp;
        const Scalar<option, errorTrack>& scalar;
    public:
        VectorExpression(const T& exp_, const Scalar<option, errorTrack>& scalar_) : exp(exp_), scalar(scalar_) {}

        typename Base::ScalarType operator[](size_t s) const { return exp[s] + scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class T1, class T2>
    class VectorExpression<Utils::ExpressionType::Sub, T1, T2>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Sub, T1, T2>> {
        static_assert(std::is_same<typename T1::VectorType, typename T2::VectorType>::value
                      , "Types of two operands of sub must be same.");
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Sub, T1, T2>>;
        T1 exp1;
        T2 exp2;
    public:
        VectorExpression(const T1& exp1_, const T2& exp2_) : exp1(exp1_), exp2(exp2_) {
            assert(exp1.getLength() == exp2.getLength());
        }

        typename Base::ScalarType operator[](size_t s) const { return exp1[s] - exp2[s]; }
        [[nodiscard]] size_t getLength() const { return exp1.getLength(); }
    };

    template<class T, ScalarOption option, bool errorTrack>
    class VectorExpression<Utils::ExpressionType::Sub, T, Scalar<option, errorTrack>>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Sub, T, Scalar<option, errorTrack>>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Sub, T, Scalar<option, errorTrack>>>;
        T exp;
        const Scalar<option, errorTrack>& scalar;
    public:
        VectorExpression(const T& exp_, const Scalar<option, errorTrack>& scalar_) : exp(exp_), scalar(scalar_) {}

        typename Base::ScalarType operator[](size_t s) const { return exp[s] - scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class T, ScalarOption option, bool errorTrack>
    class VectorExpression<Utils::ExpressionType::Mul, T, Scalar<option, errorTrack>>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Mul, T, Scalar<option, errorTrack>>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Mul, T, Scalar<option, errorTrack>>>;
        T exp;
        const Scalar<option, errorTrack>& scalar;
    public:
        VectorExpression(const T& exp_, const Scalar<option, errorTrack>& scalar_) : exp(exp_), scalar(scalar_) {}

        typename Base::ScalarType operator[](size_t s) const { return exp[s] * scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Div//////////////////////////////////////
    template<class T, ScalarOption option, bool errorTrack>
    class VectorExpression<Utils::ExpressionType::Div, T, Scalar<option, errorTrack>>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Div, T, Scalar<option, errorTrack>>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Div, T, Scalar<option, errorTrack>>>;
        T exp;
        const Scalar<option, errorTrack>& scalar;
    public:
        VectorExpression(const T& exp_, const Scalar<option, errorTrack>& scalar_) : exp(exp_), scalar(scalar_) {}

        typename Base::ScalarType operator[](size_t s) const { return exp[s] / scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class Derived>
    inline VectorExpression<Utils::ExpressionType::Minus, Derived> operator-(const VectorBase<Derived>& v) {
        return VectorExpression<Utils::ExpressionType::Minus, Derived>(v.getDerived());
    }

    template<class VectorType>
    inline VectorExpression<Utils::ExpressionType::Minus, VectorBlock<VectorType>>
    operator-(const VectorBlock<VectorType>& block) {
        return VectorExpression<Utils::ExpressionType::Minus, VectorBlock<VectorType>>(block);
    }

    template<Utils::ExpressionType type, class T1, class T2>
    inline VectorExpression<Utils::ExpressionType::Minus, VectorExpression<type, T1, T2>>
    operator-(const VectorExpression<type, T1, T2>& exp) {
        return VectorExpression<Utils::ExpressionType::Minus
                                , VectorExpression<type, T1, T2>>(exp);
    }
    //////////////////////////////////////Add//////////////////////////////////////
    template<class Derived, class OtherDerived>
    inline VectorExpression<Utils::ExpressionType::Add, Derived, OtherDerived>
            operator+(const VectorBase<Derived>& v1, const VectorBase<OtherDerived>& v2) {
        return VectorExpression<Utils::ExpressionType::Add, Derived, OtherDerived>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarType> operator+(const VectorBase<VectorType>& v, const ScalarBase<ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarType>(v.getDerived(), s.getDerived());
    }

    template<class ScalarType, class VectorType>
    inline VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarType> operator+(const ScalarBase<ScalarType>& s, const VectorBase<VectorType>& v) { return v + s; }

    template<class VectorType, class T>
    inline VectorExpression<Utils::ExpressionType::Add, VectorBlock<VectorType>, T>
    operator+(const VectorBlock<VectorType>& block, T t) {
        return VectorExpression<Utils::ExpressionType::Add, VectorBlock<VectorType>, T>(block, t);
    }

    template<Utils::ExpressionType type, class T1, class T2, ScalarOption option, bool errorTrack>
    inline VectorExpression<Utils::ExpressionType::Add, VectorExpression<type, T1, T2>, Scalar<option, errorTrack>>
    operator+(const VectorExpression<type, T1, T2>& exp, const Scalar<option, errorTrack>& s) {
        return VectorExpression<Utils::ExpressionType::Add
                                , VectorExpression<type, T1, T2>
                                , Scalar<option, errorTrack>>(exp, s);
    }

    template<Utils::ExpressionType type, class T1, class T2, class T, size_t Length, size_t MaxLength>
    inline VectorExpression<Utils::ExpressionType::Add, VectorExpression<type, T1, T2>, Vector<T, Length, MaxLength>>
    operator+(const VectorExpression<type, T1, T2>& exp, const Vector<T, Length, MaxLength>& v) {
        return VectorExpression<Utils::ExpressionType::Add
                                , VectorExpression<type, T1, T2>
                                , Vector<T, Length, MaxLength>>(exp, v);
    }

    template<Utils::ExpressionType type, class T1, class T2, class T, size_t Length, size_t MaxLength>
    inline VectorExpression<Utils::ExpressionType::Add, VectorExpression<type, T1, T2>, Vector<T, Length, MaxLength>>
    operator+(const Vector<T, Length, MaxLength>& v, const VectorExpression<type, T1, T2>& exp) {
        return exp + v;
    }

    template<Utils::ExpressionType type1, class T11, class T12, Utils::ExpressionType type2, class T21, class T22>
    inline VectorExpression<Utils::ExpressionType::Add, VectorExpression<type1, T11, T12>, VectorExpression<type2, T21, T22>>
    operator+(const VectorExpression<type1, T11, T12>& exp1, const VectorExpression<type2, T21, T22>& exp2) {
        return VectorExpression<Utils::ExpressionType::Add
                                , VectorExpression<type1, T11, T12>
                                , VectorExpression<type2, T21, T22>>(exp1, exp2);
    }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class Derived, class OtherDerived>
    inline VectorExpression<Utils::ExpressionType::Sub, Derived, OtherDerived>
            operator-(const VectorBase<Derived>& v1, const VectorBase<OtherDerived>& v2) {
        return VectorExpression<Utils::ExpressionType::Sub, Derived, OtherDerived>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    VectorExpression<Utils::ExpressionType::Sub, VectorType, ScalarType> operator+(const VectorBase<VectorType>& v, const ScalarBase<ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::Sub, VectorType, ScalarType>(v.getDerived(), s);
    }

    template<class VectorType, class T>
    inline VectorExpression<Utils::ExpressionType::Sub, VectorBlock<VectorType>, T>
    operator-(const VectorBlock<VectorType>& block, T t) {
        return VectorExpression<Utils::ExpressionType::Sub, VectorBlock<VectorType>, T>(block, t);
    }

    template<Utils::ExpressionType type, class T1, class T2, ScalarOption option, bool errorTrack>
    inline VectorExpression<Utils::ExpressionType::Sub, VectorExpression<type, T1, T2>, Scalar<option, errorTrack>>
    operator-(const VectorExpression<type, T1, T2>& exp, const Scalar<option, errorTrack>& s) {
        return VectorExpression<Utils::ExpressionType::Sub
                                , VectorExpression<type, T1, T2>
                                , Scalar<option, errorTrack>>(exp, s);
    }

    template<Utils::ExpressionType type, class T1, class T2, class T, size_t Length, size_t MaxLength>
    inline VectorExpression<Utils::ExpressionType::Sub, VectorExpression<type, T1, T2>, Vector<T, Length, MaxLength>>
    operator-(const VectorExpression<type, T1, T2>& exp, const Vector<T, Length, MaxLength>& v) {
        return VectorExpression<Utils::ExpressionType::Sub
                                , VectorExpression<type, T1, T2>
                                , Vector<T, Length, MaxLength>>(exp, v);
    }

    template<Utils::ExpressionType type, class T1, class T2, class T, size_t Length, size_t MaxLength>
    inline VectorExpression<Utils::ExpressionType::Sub, Vector<T, Length, MaxLength>, VectorExpression<type, T1, T2>>
    operator-(const Vector<T, Length, MaxLength>& v, const VectorExpression<type, T1, T2>& exp) {
        return VectorExpression<Utils::ExpressionType::Sub
                                , Vector<T, Length, MaxLength>
                                , VectorExpression<type, T1, T2>>(v, exp);
    }

    template<Utils::ExpressionType type1, class T11, class T12, Utils::ExpressionType type2, class T21, class T22>
    inline VectorExpression<Utils::ExpressionType::Sub, VectorExpression<type1, T11, T12>, VectorExpression<type2, T21, T22>>
    operator-(const VectorExpression<type1, T11, T12>& exp1, const VectorExpression<type2, T21, T22>& exp2) {
        return VectorExpression<Utils::ExpressionType::Sub
                                , VectorExpression<type1, T11, T12>
                                , VectorExpression<type2, T21, T22>>(exp1, exp2);
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class VectorType, class ScalarType>
    VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarType> operator*(const VectorBase<VectorType>& v, const ScalarBase<ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarType>(v.getDerived(), s.getDerived());
    }

    template<class ScalarType, class VectorType>
    inline VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarType> operator*(const ScalarBase<ScalarType>& s, const VectorBase<VectorType>& v) { return v * s; }

    template<class VectorType, ScalarOption option, bool errorTrack>
    inline VectorExpression<Utils::ExpressionType::Mul, VectorBlock<VectorType>, Scalar<option, errorTrack>>
    operator*(const VectorBlock<VectorType>& block, const Scalar<option, errorTrack>& s) {
        return VectorExpression<Utils::ExpressionType::Mul, VectorBlock<VectorType>, Scalar<option, errorTrack>>(block, s);
    }

    template<Utils::ExpressionType type, class T1, class T2, ScalarOption option, bool errorTrack>
    inline VectorExpression<Utils::ExpressionType::Mul, VectorExpression<type, T1, T2>, Scalar<option, errorTrack>>
    operator*(const VectorExpression<type, T1, T2>& exp, const Scalar<option, errorTrack>& s) {
        return VectorExpression<Utils::ExpressionType::Mul
                                , VectorExpression<type, T1, T2>
                                , Scalar<option, errorTrack>>(exp, s);
    }
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class ScalarType>
    VectorExpression<Utils::ExpressionType::Div, VectorType, ScalarType> operator+(const VectorBase<VectorType>& v, const ScalarBase<ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::Div, VectorType, ScalarType>(v.getDerived(), s);
    }

    template<class VectorType, ScalarOption option, bool errorTrack>
    inline VectorExpression<Utils::ExpressionType::Div, VectorBlock<VectorType>, Scalar<option, errorTrack>>
    operator/(const VectorBlock<VectorType>& block, const Scalar<option, errorTrack>& s) {
        return VectorExpression<Utils::ExpressionType::Div, VectorBlock<VectorType>, Scalar<option, errorTrack>>(block, s);
    }

    template<Utils::ExpressionType type, class T1, class T2, ScalarOption option, bool errorTrack>
    inline VectorExpression<Utils::ExpressionType::Div, VectorExpression<type, T1, T2>, Scalar<option, errorTrack>>
    operator/(const VectorExpression<type, T1, T2>& exp, const Scalar<option, errorTrack>& s) {
        return VectorExpression<Utils::ExpressionType::Div
                                , VectorExpression<type, T1, T2>
                                , Scalar<option, errorTrack>>(exp, s);
    }
}
