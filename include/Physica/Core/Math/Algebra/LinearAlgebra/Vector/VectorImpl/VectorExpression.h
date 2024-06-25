/*
 * Copyright 2020-2024 WeiBo He.
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
#include "Physica/Core/MultiPrecision/ScalarImpl/ExpressionType.h"
#include "Physica/Utils/Template/CRTPBase.h"

namespace Physica::Core {
    //Forward declaration
    using Utils::Dynamic;
    /**
     * \class VectorExpression implements template expression for vectors, which will reduce temporary objects.
     * 
     * Operations defined as \tparam T1 \tparam Type \tparam T2. e.g. vector + scalar, expression * expression
     */
    template<ExpressionType Type, class T1, class T2 = T1>
    class VectorExpression;

    namespace Internal {
        template<class T>
        class Traits;

        template<ExpressionType Type, class Expr1, class Expr2>
        class Traits<VectorExpression<Type, Expr1, Expr2>> {
            static_assert(Expr1::SizeAtCompile == Dynamic || Expr2::SizeAtCompile == Dynamic || (Expr1::SizeAtCompile == Expr2::SizeAtCompile)
                         , "[Error]: Vector dimentions do not match");
            using ScalarType1 = typename Expr1::ScalarType;
            using RealType = typename ScalarType1::RealType;
            using BinaryScalarType = typename BinaryScalarOpReturnType<ScalarType1, typename Expr2::ScalarType>::Type;
            constexpr static bool FastAssign1 = Traits<Expr1>::FastAssign;
            constexpr static bool FastAssign2 = Traits<Expr2>::FastAssign;
            constexpr static bool IsAddOrSub = Type == ExpressionType::Add || Type == ExpressionType::Sub;
        public:
            using ScalarType = typename std::conditional<Type == ExpressionType::Abs, RealType, BinaryScalarType>::type;
            constexpr static size_t SizeAtCompile = Expr1::SizeAtCompile > Expr2::SizeAtCompile ? Expr1::SizeAtCompile : Expr2::SizeAtCompile;
            constexpr static size_t MaxSizeAtCompile = Expr1::MaxSizeAtCompile > Expr2::MaxSizeAtCompile ? Expr1::MaxSizeAtCompile : Expr2::MaxSizeAtCompile;

            using PacketType = typename BestPacket<ScalarType, SizeAtCompile>::Type;
            constexpr static bool FastAssign = IsAddOrSub && (FastAssign1 || FastAssign2);
        };

        template<ExpressionType Type, class Expr, class AnyScalar>
        class Traits<VectorExpression<Type, Expr, ScalarBase<AnyScalar>>> {
            static_assert(is_scalar<AnyScalar>::value, "[Error]: This is not a scalar type");
        public:
            using ScalarType = typename BinaryScalarOpReturnType<typename Expr::ScalarType, AnyScalar>::Type;
            constexpr static size_t SizeAtCompile = Expr::SizeAtCompile;
            constexpr static size_t MaxSizeAtCompile = Expr::MaxSizeAtCompile;

            using PacketType = typename BestPacket<ScalarType, SizeAtCompile>::Type;
            constexpr static bool FastAssign = Traits<Expr>::FastAssign;
        };
    }
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class VectorType>
    class VectorExpression<ExpressionType::Minus, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Minus, VectorType>> {
        using This = VectorExpression<ExpressionType::Minus, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& exp;
    public:
        explicit VectorExpression(const RValueVector<VectorType>& exp_) : exp(exp_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return -exp.calc(s); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const { return -exp.template packet<PacketType>(index); }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const { return -exp.template packetPartial<PacketType>(index, count); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class AnyScalar>
    class VectorExpression<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpression<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>> {
        using This = VectorExpression<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>;
        using Base = RValueVector<This>;
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return exp.calc(s) / scalar; }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp.template packet<PacketType>(index) * PacketType(reciprocal(scalar));
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const {
            return exp.template packetPartial<PacketType>(index, count) * PacketType(reciprocal(scalar));
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType1, class VectorType2>
    class VectorExpression<ExpressionType::Div, VectorType1, VectorType2>
            : public RValueVector<VectorExpression<ExpressionType::Div, VectorType1, VectorType2>> {
        using This = VectorExpression<ExpressionType::Div, VectorType1, VectorType2>;
        using Base = RValueVector<This>;
        const VectorType1& v1;
        const VectorType2& v2;
    public:
        VectorExpression(const RValueVector<VectorType1>& v1_, const RValueVector<VectorType2>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1.getLength() == v2.getLength());
        }
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return v1.calc(s) / v2.calc(s); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return v1.template packet<PacketType>(index) / v2.template packet<PacketType>(index);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const {
            return v1.template packetPartial<PacketType>(index, count) / v2.template packetPartial<PacketType>(index, count);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v1.getLength(); }
    };
    //////////////////////////////////////Compare//////////////////////////////////////
    template<class VectorType1, class VectorType2>
    class VectorExpression<ExpressionType::More, VectorType1, VectorType2>
            : public RValueVector<VectorExpression<ExpressionType::More, VectorType1, VectorType2>> {
        using This = VectorExpression<ExpressionType::More, VectorType1, VectorType2>;
        using Base = RValueVector<This>;
        const VectorType1& v1;
        const VectorType2& v2;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpression(const RValueVector<VectorType1>& v1_, const RValueVector<VectorType2>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1.getLength() == v2.getLength());
        }
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(v1.calc(s) > v2.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return v1.template packet<PacketType>(index) > v2.template packet<PacketType>(index);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const {
            return v1.template packetPartial<PacketType>(index, count) > v2.template packetPartial<PacketType>(index, count);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v1.getLength(); }
    };

    template<class VectorType, class AnyScalar>
    class VectorExpression<ExpressionType::More, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpression<ExpressionType::More, VectorType, ScalarBase<AnyScalar>>> {
        using This = VectorExpression<ExpressionType::More, VectorType, ScalarBase<AnyScalar>>;
        using Base = RValueVector<This>;
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(exp.calc(s) > scalar); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp.template packet<PacketType>(index) > PacketType(scalar);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const {
            return exp.template packetPartial<PacketType>(index, count) > PacketType(scalar);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType1, class VectorType2>
    class VectorExpression<ExpressionType::MoreEq, VectorType1, VectorType2>
            : public RValueVector<VectorExpression<ExpressionType::MoreEq, VectorType1, VectorType2>> {
        using This = VectorExpression<ExpressionType::MoreEq, VectorType1, VectorType2>;
        using Base = RValueVector<This>;
        const VectorType1& v1;
        const VectorType2& v2;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpression(const RValueVector<VectorType1>& v1_, const RValueVector<VectorType2>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1.getLength() == v2.getLength());
        }
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(v1.calc(s) >= v2.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return v1.template packet<PacketType>(index) >= v2.template packet<PacketType>(index);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const {
            return v1.template packetPartial<PacketType>(index, count) >= v2.template packetPartial<PacketType>(index, count);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v1.getLength(); }
    };

    template<class VectorType, class AnyScalar>
    class VectorExpression<ExpressionType::MoreEq, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpression<ExpressionType::MoreEq, VectorType, ScalarBase<AnyScalar>>> {
        using This = VectorExpression<ExpressionType::MoreEq, VectorType, ScalarBase<AnyScalar>>;
        using Base = RValueVector<This>;
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(exp.calc(s) >= scalar); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp.template packet<PacketType>(index) >= PacketType(scalar);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const {
            return exp.template packetPartial<PacketType>(index, count) >= PacketType(scalar);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    class VectorExpression<ExpressionType::Reciprocal, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Reciprocal, VectorType>> {
        using This = VectorExpression<ExpressionType::Reciprocal, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& exp;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_) : exp(exp_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t index) const { return reciprocal(exp.calc(index)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const { return PacketType(1) / exp.template packet<PacketType>(index); }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const { return PacketType(1) / exp.template packetPartial<PacketType>(index, count); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<ExpressionType::Sqrt, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Sqrt, VectorType>> {
        using This = VectorExpression<ExpressionType::Sqrt, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& exp;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_) : exp(exp_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return sqrt(exp.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const { return sqrt(exp.template packet<PacketType>(index)); }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const { return sqrt(exp.template packetPartial<PacketType>(index, count)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<ExpressionType::Cbrt, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Cbrt, VectorType>> {
        using This = VectorExpression<ExpressionType::Cbrt, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& exp;
    public:
        VectorExpression(const RValueVector<VectorType>& exp_) : exp(exp_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return cbrt(exp.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            PacketType result = exp.template packet<PacketType>(index);
            for (size_t i = 0; i < static_cast<size_t>(PacketType::size()); ++i)
                result.insert(i, cbrt(result[i]));
            return result;
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const {
            PacketType result = exp.template packetPartial<PacketType>(index, count);
            for (size_t i = 0; i < count; ++i)
                result.insert(i, cbrt(result[i]));
            return result;
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<ExpressionType::Abs, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Abs, VectorType>> {
        using This = VectorExpression<ExpressionType::Abs, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return abs(v.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const { return abs(v.template packet<PacketType>(index)); }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const { return abs(v.template packetPartial<PacketType>(index, count)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<ExpressionType::Square, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Square, VectorType>> {
        using This = VectorExpression<ExpressionType::Square, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return square(v.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const { return square(v.template packet<PacketType>(index)); }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const { return square(v.template packetPartial<PacketType>(index, count)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<ExpressionType::Ln, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Ln, VectorType>> {
        using This = VectorExpression<ExpressionType::Ln, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return ln(v.calc(s)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<ExpressionType::Exp, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Exp, VectorType>> {
        using This = VectorExpression<ExpressionType::Exp, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return exp(v.calc(s)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<ExpressionType::Pow, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Pow, VectorType>> {
        using This = VectorExpression<ExpressionType::Pow, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
        const ScalarType& s;
    public:
        VectorExpression(const RValueVector<VectorType>& v_, const ScalarBase<ScalarType>& s_)
                : v(v_.getDerived()), s(s_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return pow(v.calc(i), s); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<ExpressionType::Sin, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Sin, VectorType>> {
        using This = VectorExpression<ExpressionType::Sin, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return sin(v.calc(i)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<ExpressionType::Cos, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Cos, VectorType>> {
        using This = VectorExpression<ExpressionType::Cos, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return cos(v.calc(i)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<ExpressionType::Relu, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Relu, VectorType>> {
        using This = VectorExpression<ExpressionType::Relu, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return relu(v.calc(i)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpression<ExpressionType::Unit, VectorType>
            : public RValueVector<VectorExpression<ExpressionType::Unit, VectorType>> {
        using This = VectorExpression<ExpressionType::Unit, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
    public:
        VectorExpression(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return v.calc(i).unit(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    template<class Derived>
    [[nodiscard]] inline VectorExpression<ExpressionType::Minus, Derived> operator-(const RValueVector<Derived>& v) noexcept {
        return VectorExpression<ExpressionType::Minus, Derived>(v.getDerived());
    }
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class ScalarType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Div, VectorType, ScalarBase<ScalarType>>
    operator/(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) noexcept {
        return VectorExpression<ExpressionType::Div, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline VectorExpression<ExpressionType::Div, VectorType1, VectorType2>
    divide(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpression<ExpressionType::Div, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }
    //////////////////////////////////////Compare//////////////////////////////////////
    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline VectorExpression<ExpressionType::More, VectorType1, VectorType2>
    operator>(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpression<ExpressionType::More, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline VectorExpression<ExpressionType::MoreEq, VectorType1, VectorType2>
    operator>=(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpression<ExpressionType::MoreEq, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    [[nodiscard]] inline VectorExpression<ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>
    operator>=(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) noexcept {
        return VectorExpression<ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline VectorExpression<ExpressionType::More, VectorType2, VectorType1>
    operator<(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpression<ExpressionType::More, VectorType2, VectorType1>(v2.getDerived(), v1.getDerived());
    }

    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline VectorExpression<ExpressionType::MoreEq, VectorType2, VectorType1>
    operator<=(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpression<ExpressionType::MoreEq, VectorType2, VectorType1>(v2.getDerived(), v1.getDerived());
    }

    template<class VectorType, class ScalarType>
    [[nodiscard]] inline VectorExpression<ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>
    operator<=(const ScalarBase<ScalarType>& s, const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Reciprocal, VectorType>
    reciprocal(const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::Reciprocal, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Sqrt, VectorType> sqrt(const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::Sqrt, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Cbrt, VectorType> cbrt(const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::Cbrt, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Abs, VectorType> abs(const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::Abs, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Relu, VectorType> relu(
            const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::Relu, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Unit, VectorType> unit(const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::Unit, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Square, VectorType> square(const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::Square, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Ln, VectorType> ln(const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::Ln, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Exp, VectorType> exp(const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::Exp, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Pow, VectorType> pow(
            const RValueVector<VectorType>& v,
            const ScalarBase<typename VectorType::ScalarType>& s) noexcept {
        return VectorExpression<ExpressionType::Pow, VectorType>(v, s);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Sin, VectorType> sin(
            const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::Sin, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Cos, VectorType> cos(
            const RValueVector<VectorType>& v) noexcept {
        return VectorExpression<ExpressionType::Cos, VectorType>(v);
    }
}

#include "VectorExprImpl/VectorAdd.h"
#include "VectorExprImpl/VectorSub.h"
#include "VectorExprImpl/VectorMul.h"
