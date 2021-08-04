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
        public:
            using ScalarType = typename BinaryScalarOpReturnType<typename Exp1::ScalarType, typename Exp2::ScalarType>::Type;
        };
        /**
         * This class implements assignTo() for all \class VectorExpression
         */
        template<class Derived>
        class VectorExpressionHelper : public RValueVector<Derived> {
        public:
            using ScalarType = typename Traits<Derived>::ScalarType;
        private:
            using Base = RValueVector<Derived>;
        public:
            template<class OtherDerived>
            void assignTo(LValueVector<OtherDerived>& v) const {
                assert(v.getLength() == getLength());
                for (size_t i = 0; i < getLength(); ++i)
                    v[i] = calc(i);
            }

            [[nodiscard]] ScalarType calc(size_t index) const { return Base::getDerived().calc(index); }
            /* Getters */
            [[nodiscard]] size_t getLength() const noexcept { return Base::getDerived().getLength(); }
        };
    }
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class VectorType>
    class VectorExpression<Utils::ExpressionType::Minus, VectorType>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Minus, VectorType>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Minus, VectorType>>;
        const VectorType& exp;
    public:
        explicit VectorExpression(const RValueVector<VectorType>& exp_) : exp(exp_.getDerived()) {}

        typename Base::ScalarType calc(size_t s) const { return -exp.calc(s); }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Add//////////////////////////////////////
    template<class VectorType1, class VectorType2>
    class VectorExpression<Utils::ExpressionType::Add, VectorType1, VectorType2>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Add, VectorType1, VectorType2>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Add, VectorType1, VectorType2>>;
        const VectorType1& exp1;
        const VectorType2& exp2;
    public:
        VectorExpression(const RValueVector<VectorType1>& exp1_, const RValueVector<VectorType2>& exp2_)
                : exp1(exp1_.getDerived()), exp2(exp2_.getDerived()) {
            assert(exp1.getLength() == exp2.getLength());
        }

        typename Base::ScalarType calc(size_t s) const { return exp1.calc(s) + exp2.calc(s); }
        [[nodiscard]] size_t getLength() const { return exp1.getLength(); }
    };

    template<class VectorType, ScalarOption option, bool errorTrack>
    class VectorExpression<Utils::ExpressionType::Add, VectorType, Scalar<option, errorTrack>>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Add, VectorType, Scalar<option, errorTrack>>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Add, VectorType, Scalar<option, errorTrack>>>;
        const VectorType& exp;
        const Scalar<option, errorTrack>& scalar;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const Scalar<option, errorTrack>& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}

        typename Base::ScalarType calc(size_t s) const { return exp.calc(s) + scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class VectorType1, class VectorType2>
    class VectorExpression<Utils::ExpressionType::Sub, VectorType1, VectorType2>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Sub, VectorType1, VectorType2>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Sub, VectorType1, VectorType2>>;
        const VectorType1& exp1;
        const VectorType2& exp2;
    public:
        VectorExpression(const RValueVector<VectorType1>& exp1_, const RValueVector<VectorType2>& exp2_)
                : exp1(exp1_.getDerived()), exp2(exp2_.getDerived()) {
            assert(exp1.getLength() == exp2.getLength());
        }

        typename Base::ScalarType calc(size_t s) const { return exp1.calc(s) - exp2.calc(s); }
        [[nodiscard]] size_t getLength() const { return exp1.getLength(); }
    };

    template<class VectorType, ScalarOption option, bool errorTrack>
    class VectorExpression<Utils::ExpressionType::Sub, VectorType, Scalar<option, errorTrack>>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Sub, VectorType, Scalar<option, errorTrack>>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Sub, VectorType, Scalar<option, errorTrack>>>;
        const VectorType& exp;
        const Scalar<option, errorTrack>& scalar;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const Scalar<option, errorTrack>& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}

        typename Base::ScalarType calc(size_t s) const { return exp.calc(s) - scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class VectorType, class ScalarType>
    class VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarType>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarType>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarType>>;
        const VectorType& exp;
        const ScalarType& scalar;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const ScalarBase<ScalarType>& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_.getDerived()) {}

        typename Base::ScalarType calc(size_t s) const { return exp.calc(s) * scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class ScalarType>
    class VectorExpression<Utils::ExpressionType::Div, VectorType, ScalarType>
            : public Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Div, VectorType, ScalarType>> {
        using Base = Internal::VectorExpressionHelper<VectorExpression<Utils::ExpressionType::Div, VectorType, ScalarType>>;
        const VectorType& exp;
        const ScalarType& scalar;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const ScalarBase<ScalarType>& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_.getDerived()) {}

        typename Base::ScalarType calc(size_t s) const { return exp.calc(s) / scalar; }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class Derived>
    inline VectorExpression<Utils::ExpressionType::Minus, Derived> operator-(const RValueVector<Derived>& v) {
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
            operator+(const RValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        return VectorExpression<Utils::ExpressionType::Add, Derived, OtherDerived>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarType> operator+(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarType>(v.getDerived(), s.getDerived());
    }

    template<class ScalarType, class VectorType>
    inline VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarType> operator+(const ScalarBase<ScalarType>& s, const RValueVector<VectorType>& v) { return v + s; }

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
            operator-(const RValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        return VectorExpression<Utils::ExpressionType::Sub, Derived, OtherDerived>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    VectorExpression<Utils::ExpressionType::Sub, VectorType, ScalarType> operator+(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) {
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
    VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarType> operator*(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarType>(v.getDerived(), s.getDerived());
    }

    template<class ScalarType, class VectorType>
    inline VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarType> operator*(const ScalarBase<ScalarType>& s, const RValueVector<VectorType>& v) { return v * s; }

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
    VectorExpression<Utils::ExpressionType::Div, VectorType, ScalarType> operator+(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) {
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
