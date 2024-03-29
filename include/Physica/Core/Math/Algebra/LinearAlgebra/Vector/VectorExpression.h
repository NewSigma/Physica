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

        template<Utils::ExpressionType type, class Exp1, class Exp2>
        class Traits<VectorExpression<type, Exp1, Exp2>> {
            static_assert(Exp1::SizeAtCompile == Exp2::SizeAtCompile, "[Possible optimization]: Binary operation between vectors should have same size at compile");
            static_assert(Exp1::MaxSizeAtCompile == Exp2::MaxSizeAtCompile, "[Possible optimization]: Binary operation between vectors should have same max size at compile");
        public:
            using ScalarType = typename BinaryScalarOpReturnType<typename Exp1::ScalarType, typename Exp2::ScalarType>::Type;
            constexpr static size_t SizeAtCompile = Exp1::SizeAtCompile;
            constexpr static size_t MaxSizeAtCompile = Exp1::MaxSizeAtCompile;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };

        template<Utils::ExpressionType type, class Exp, class AnyScalar>
        class Traits<VectorExpression<type, Exp, ScalarBase<AnyScalar>>> {
        public:
            using ScalarType = typename BinaryScalarOpReturnType<typename Exp::ScalarType, AnyScalar>::Type;
            constexpr static size_t SizeAtCompile = Exp::SizeAtCompile;
            constexpr static size_t MaxSizeAtCompile = Exp::MaxSizeAtCompile;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };
    }
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class VectorType>
    class VectorExpression<Utils::ExpressionType::Minus, VectorType>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Minus, VectorType>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Minus, VectorType>>;
        const VectorType& exp;
    public:
        explicit VectorExpression(const RValueVector<VectorType>& exp_) : exp(exp_.getDerived()) {}

        typename Base::ScalarType calc(size_t s) const { return -exp.calc(s); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const { return -exp.template packet<PacketType>(index); }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const { return -exp.template packetPartial<PacketType>(index); }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Add//////////////////////////////////////
    template<class VectorType1, class VectorType2>
    class VectorExpression<Utils::ExpressionType::Add, VectorType1, VectorType2>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Add, VectorType1, VectorType2>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Add, VectorType1, VectorType2>>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType1& exp1;
        const VectorType2& exp2;
    public:
        VectorExpression(const RValueVector<VectorType1>& exp1_, const RValueVector<VectorType2>& exp2_)
                : exp1(exp1_.getDerived()), exp2(exp2_.getDerived()) {
            assert(exp1.getLength() == exp2.getLength());
        }

        ScalarType calc(size_t s) const { return ScalarType(exp1.calc(s)) + ScalarType(exp2.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp1.template packet<PacketType>(index) + exp2.template packet<PacketType>(index);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return exp1.template packetPartial<PacketType>(index) + exp2.template packetPartial<PacketType>(index);
        }
        [[nodiscard]] size_t getLength() const { return exp1.getLength(); }
    };

    template<class VectorType, class AnyScalar>
    class VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarBase<AnyScalar>>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarBase<AnyScalar>>>;
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}

        ScalarType calc(size_t s) const { return ScalarType(exp.calc(s)) + ScalarType(scalar); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp.template packet<PacketType>(index) + PacketType(scalar.getTrivial());
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return exp.template packetPartial<PacketType>(index) + PacketType(scalar.getTrivial());
        }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class VectorType1, class VectorType2>
    class VectorExpression<Utils::ExpressionType::Sub, VectorType1, VectorType2>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Sub, VectorType1, VectorType2>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Sub, VectorType1, VectorType2>>;
        const VectorType1& exp1;
        const VectorType2& exp2;
    public:
        VectorExpression(const RValueVector<VectorType1>& exp1_, const RValueVector<VectorType2>& exp2_)
                : exp1(exp1_.getDerived()), exp2(exp2_.getDerived()) {
            assert(exp1.getLength() == exp2.getLength());
        }

        typename Base::ScalarType calc(size_t s) const { return exp1.calc(s) - exp2.calc(s); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp1.template packet<PacketType>(index) - exp2.template packet<PacketType>(index);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return exp1.template packetPartial<PacketType>(index) - exp2.template packetPartial<PacketType>(index);
        }
        [[nodiscard]] size_t getLength() const { return exp1.getLength(); }
    };

    template<class VectorType, class AnyScalar>
    class VectorExpression<Utils::ExpressionType::Sub, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Sub, VectorType, ScalarBase<AnyScalar>>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Sub, VectorType, ScalarBase<AnyScalar>>>;
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}

        typename Base::ScalarType calc(size_t s) const { return exp.calc(s) - scalar; }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp.template packet<PacketType>(index) - PacketType(scalar.getTrivial());
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return exp.template packetPartial<PacketType>(index) - PacketType(scalar.getTrivial());
        }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class VectorType, class AnyScalar>
    class VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarBase<AnyScalar>>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarBase<AnyScalar>>>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}

        ScalarType calc(size_t s) const { return ScalarType(exp.calc(s)) * ScalarType(scalar); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp.template packet<PacketType>(index) * PacketType(scalar.getTrivial());
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return exp.template packetPartial<PacketType>(index) * PacketType(scalar.getTrivial());
        }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType1, class VectorType2>
    class VectorExpression<Utils::ExpressionType::Mul, VectorType1, VectorType2>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Mul, VectorType1, VectorType2>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Mul, VectorType1, VectorType2>>;
        const VectorType1& v1;
        const VectorType2& v2;
    public:
        VectorExpression(const RValueVector<VectorType1>& v1_, const RValueVector<VectorType2>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1.getLength() == v2.getLength());
        }

        typename Base::ScalarType calc(size_t s) const { return v1.calc(s) * v2.calc(s); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return v1.template packet<PacketType>(index) * v2.template packet<PacketType>(index);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return v1.template packetPartial<PacketType>(index) * v2.template packetPartial<PacketType>(index);
        }
        [[nodiscard]] size_t getLength() const { return v1.getLength(); }
    };
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class AnyScalar>
    class VectorExpression<Utils::ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>>;
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}

        typename Base::ScalarType calc(size_t s) const { return exp.calc(s) / scalar; }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp.template packet<PacketType>(index) * PacketType(reciprocal(scalar).getTrivial());
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return exp.template packetPartial<PacketType>(index) * PacketType(reciprocal(scalar).getTrivial());
        }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType1, class VectorType2>
    class VectorExpression<Utils::ExpressionType::Div, VectorType1, VectorType2>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Div, VectorType1, VectorType2>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Div, VectorType1, VectorType2>>;
        const VectorType1& v1;
        const VectorType2& v2;
    public:
        VectorExpression(const RValueVector<VectorType1>& v1_, const RValueVector<VectorType2>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1.getLength() == v2.getLength());
        }

        typename Base::ScalarType calc(size_t s) const { return v1.calc(s) / v2.calc(s); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return v1.template packet<PacketType>(index) / v2.template packet<PacketType>(index);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return v1.template packetPartial<PacketType>(index) / v2.template packetPartial<PacketType>(index);
        }
        [[nodiscard]] size_t getLength() const { return v1.getLength(); }
    };
    //////////////////////////////////////Compare//////////////////////////////////////
    template<class VectorType1, class VectorType2>
    class VectorExpression<Utils::ExpressionType::More, VectorType1, VectorType2>
            : public RValueVector<VectorExpression<Utils::ExpressionType::More, VectorType1, VectorType2>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::More, VectorType1, VectorType2>>;
        const VectorType1& v1;
        const VectorType2& v2;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpression(const RValueVector<VectorType1>& v1_, const RValueVector<VectorType2>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1.getLength() == v2.getLength());
        }

        ScalarType calc(size_t s) const { return ScalarType(v1.calc(s) > v2.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return v1.template packet<PacketType>(index) > v2.template packet<PacketType>(index);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return v1.template packetPartial<PacketType>(index) > v2.template packetPartial<PacketType>(index);
        }
        [[nodiscard]] size_t getLength() const { return v1.getLength(); }
    };

    template<class VectorType, class AnyScalar>
    class VectorExpression<Utils::ExpressionType::More, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpression<Utils::ExpressionType::More, VectorType, ScalarBase<AnyScalar>>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::More, VectorType, ScalarBase<AnyScalar>>>;
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}

        ScalarType calc(size_t s) const { return ScalarType(exp.calc(s) > scalar); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp.template packet<PacketType>(index) > PacketType(scalar.getTrivial());
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return exp.template packetPartial<PacketType>(index) > PacketType(scalar.getTrivial());
        }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType1, class VectorType2>
    class VectorExpression<Utils::ExpressionType::MoreEq, VectorType1, VectorType2>
            : public RValueVector<VectorExpression<Utils::ExpressionType::MoreEq, VectorType1, VectorType2>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::MoreEq, VectorType1, VectorType2>>;
        const VectorType1& v1;
        const VectorType2& v2;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpression(const RValueVector<VectorType1>& v1_, const RValueVector<VectorType2>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1.getLength() == v2.getLength());
        }

        ScalarType calc(size_t s) const { return ScalarType(v1.calc(s) >= v2.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return v1.template packet<PacketType>(index) >= v2.template packet<PacketType>(index);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return v1.template packetPartial<PacketType>(index) >= v2.template packetPartial<PacketType>(index);
        }
        [[nodiscard]] size_t getLength() const { return v1.getLength(); }
    };

    template<class VectorType, class AnyScalar>
    class VectorExpression<Utils::ExpressionType::MoreEq, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpression<Utils::ExpressionType::MoreEq, VectorType, ScalarBase<AnyScalar>>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::MoreEq, VectorType, ScalarBase<AnyScalar>>>;
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}

        ScalarType calc(size_t s) const { return ScalarType(exp.calc(s) >= scalar); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp.template packet<PacketType>(index) >= PacketType(scalar.getTrivial());
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const {
            return exp.template packetPartial<PacketType>(index) >= PacketType(scalar.getTrivial());
        }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    class VectorExpression<Utils::ExpressionType::Reciprocal, VectorType>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Reciprocal, VectorType>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Reciprocal, VectorType>>;
        const VectorType& exp;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_) : exp(exp_.getDerived()) {}

        typename Base::ScalarType calc(size_t s) const { return reciprocal(exp.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const { return PacketType(1) / exp.template packet<PacketType>(index); }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const { return PacketType(1) / exp.template packetPartial<PacketType>(index); }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<Utils::ExpressionType::Sqrt, VectorType>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Sqrt, VectorType>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Sqrt, VectorType>>;
        const VectorType& exp;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_) : exp(exp_.getDerived()) {}

        typename Base::ScalarType calc(size_t s) const { return sqrt(exp.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const { return sqrt(exp.template packet<PacketType>(index)); }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const { return sqrt(exp.template packetPartial<PacketType>(index)); }
        [[nodiscard]] size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<Utils::ExpressionType::Abs, VectorType>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Abs, VectorType>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Abs, VectorType>>;
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}

        typename Base::ScalarType calc(size_t s) const { return abs(v.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const { return abs(v.template packet<PacketType>(index)); }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const { return abs(v.template packetPartial<PacketType>(index)); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<Utils::ExpressionType::Square, VectorType>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Square, VectorType>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Square, VectorType>>;
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}

        typename Base::ScalarType calc(size_t s) const { return square(v.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const { return square(v.template packet<PacketType>(index)); }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const { return square(v.template packetPartial<PacketType>(index)); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<Utils::ExpressionType::Ln, VectorType>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Ln, VectorType>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Ln, VectorType>>;
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}

        typename Base::ScalarType calc(size_t s) const { return ln(v.calc(s)); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<Utils::ExpressionType::Exp, VectorType>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Exp, VectorType>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Exp, VectorType>>;
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}

        typename Base::ScalarType calc(size_t s) const { return exp(v.calc(s)); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<Utils::ExpressionType::Pow, VectorType>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Pow, VectorType>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Pow, VectorType>>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
        const ScalarType& s;
    public:
        VectorExpression(const RValueVector<VectorType>& v_, const ScalarBase<ScalarType>& s_)
                : v(v_.getDerived()), s(s_.getDerived()) {}

        ScalarType calc(size_t i) const { return pow(v.calc(i), s); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<Utils::ExpressionType::Sin, VectorType>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Sin, VectorType>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Sin, VectorType>>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}

        ScalarType calc(size_t i) const { return sin(v.calc(i)); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<Utils::ExpressionType::Cos, VectorType>
            : public RValueVector<VectorExpression<Utils::ExpressionType::Cos, VectorType>> {
        using Base = RValueVector<VectorExpression<Utils::ExpressionType::Cos, VectorType>>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}

        ScalarType calc(size_t i) const { return cos(v.calc(i)); }
        [[nodiscard]] size_t getLength() const { return v.getLength(); }
    };


    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class Derived>
    inline VectorExpression<Utils::ExpressionType::Minus, Derived> operator-(const RValueVector<Derived>& v) {
        return VectorExpression<Utils::ExpressionType::Minus, Derived>(v.getDerived());
    }
    //////////////////////////////////////Add//////////////////////////////////////
    template<class Derived, class OtherDerived>
    inline VectorExpression<Utils::ExpressionType::Add, Derived, OtherDerived>
            operator+(const RValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        return VectorExpression<Utils::ExpressionType::Add, Derived, OtherDerived>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    inline VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarBase<ScalarType>> operator+(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class ScalarType, class VectorType>
    inline VectorExpression<Utils::ExpressionType::Add, VectorType, ScalarBase<ScalarType>> operator+(const ScalarBase<ScalarType>& s, const RValueVector<VectorType>& v) { return v + s; }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class Derived, class OtherDerived>
    inline VectorExpression<Utils::ExpressionType::Sub, Derived, OtherDerived>
            operator-(const RValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        return VectorExpression<Utils::ExpressionType::Sub, Derived, OtherDerived>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    inline VectorExpression<Utils::ExpressionType::Sub, VectorType, ScalarBase<ScalarType>> operator-(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::Sub, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class VectorType, class ScalarType>
    inline VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarBase<ScalarType>> operator*(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class ScalarType, class VectorType>
    inline VectorExpression<Utils::ExpressionType::Mul, VectorType, ScalarBase<ScalarType>> operator*(const ScalarBase<ScalarType>& s, const RValueVector<VectorType>& v) { return v * s; }
    
    template<class VectorType1, class VectorType2>
    inline VectorExpression<Utils::ExpressionType::Mul, VectorType1, VectorType2> hadamard(
            const RValueVector<VectorType1>& v1,
            const RValueVector<VectorType2>& v2) {
        return VectorExpression<Utils::ExpressionType::Mul, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class ScalarType>
    inline VectorExpression<Utils::ExpressionType::Div, VectorType, ScalarBase<ScalarType>> operator/(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::Div, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class VectorType1, class VectorType2>
    inline VectorExpression<Utils::ExpressionType::Div, VectorType1, VectorType2> divide(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) {
        return VectorExpression<Utils::ExpressionType::Div, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }
    //////////////////////////////////////Compare//////////////////////////////////////
    template<class VectorType1, class VectorType2>
    inline VectorExpression<Utils::ExpressionType::More, VectorType1, VectorType2>
    operator>(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) {
        return VectorExpression<Utils::ExpressionType::More, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType1, class VectorType2>
    inline VectorExpression<Utils::ExpressionType::MoreEq, VectorType1, VectorType2>
    operator>=(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) {
        return VectorExpression<Utils::ExpressionType::MoreEq, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    inline VectorExpression<Utils::ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>
    operator>=(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class VectorType1, class VectorType2>
    inline VectorExpression<Utils::ExpressionType::More, VectorType2, VectorType1>
    operator<(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) {
        return VectorExpression<Utils::ExpressionType::More, VectorType2, VectorType1>(v2.getDerived(), v1.getDerived());
    }

    template<class VectorType1, class VectorType2>
    inline VectorExpression<Utils::ExpressionType::MoreEq, VectorType2, VectorType1>
    operator<=(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) {
        return VectorExpression<Utils::ExpressionType::MoreEq, VectorType2, VectorType1>(v2.getDerived(), v1.getDerived());
    }

    template<class VectorType, class ScalarType>
    inline VectorExpression<Utils::ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>
    operator<=(const ScalarBase<ScalarType>& s, const RValueVector<VectorType>& v) {
        return VectorExpression<Utils::ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    inline VectorExpression<Utils::ExpressionType::Reciprocal, VectorType> reciprocal(const RValueVector<VectorType>& v) {
        return VectorExpression<Utils::ExpressionType::Reciprocal, VectorType>(v);
    }

    template<class VectorType>
    inline VectorExpression<Utils::ExpressionType::Sqrt, VectorType> sqrt(const RValueVector<VectorType>& v) {
        return VectorExpression<Utils::ExpressionType::Sqrt, VectorType>(v);
    }

    template<class VectorType>
    inline VectorExpression<Utils::ExpressionType::Abs, VectorType> abs(const RValueVector<VectorType>& v) {
        return VectorExpression<Utils::ExpressionType::Abs, VectorType>(v);
    }

    template<class VectorType>
    inline VectorExpression<Utils::ExpressionType::Square, VectorType> square(const RValueVector<VectorType>& v) {
        return VectorExpression<Utils::ExpressionType::Square, VectorType>(v);
    }

    template<class VectorType>
    inline VectorExpression<Utils::ExpressionType::Ln, VectorType> ln(const RValueVector<VectorType>& v) {
        return VectorExpression<Utils::ExpressionType::Ln, VectorType>(v);
    }

    template<class VectorType>
    inline VectorExpression<Utils::ExpressionType::Exp, VectorType> exp(const RValueVector<VectorType>& v) {
        return VectorExpression<Utils::ExpressionType::Exp, VectorType>(v);
    }

    template<class VectorType>
    inline VectorExpression<Utils::ExpressionType::Pow, VectorType> pow(
            const RValueVector<VectorType>& v,
            const ScalarBase<typename VectorType::ScalarType>& s) {
        return VectorExpression<Utils::ExpressionType::Pow, VectorType>(v, s);
    }

    template<class VectorType>
    inline VectorExpression<Utils::ExpressionType::Sin, VectorType> sin(
            const RValueVector<VectorType>& v) {
        return VectorExpression<Utils::ExpressionType::Sin, VectorType>(v);
    }

    template<class VectorType>
    inline VectorExpression<Utils::ExpressionType::Cos, VectorType> cos(
            const RValueVector<VectorType>& v) {
        return VectorExpression<Utils::ExpressionType::Cos, VectorType>(v);
    }
}
