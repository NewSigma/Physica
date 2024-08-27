/*
 * Copyright 2020-2024 Weibo He.
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
#include <Physica/Core/MultiPrecision/Scalar.h>
#include <Physica/Core/MultiPrecision/ScalarImpl/ExpressionType.h>
#include <Physica/Utils/Template/CRTPBase.h>

namespace Physica::Core {
    /**
     * \class VectorExpr implements template expression for vectors, which will reduce temporary objects.
     * 
     * Operations defined as \tparam T1 \tparam Type \tparam T2. e.g. vector + scalar, expression * expression
     */
    template<ExpressionType Type, class T1, class T2 = T1> class VectorExpr;
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class AnyScalar>
    class VectorExpr<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpr<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>> {
        using This = VectorExpr<ExpressionType::Div, VectorType, ScalarBase<AnyScalar>>;
        using Base = RValueVector<This>;
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        VectorExpr(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return exp.calc(s) / scalar; }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return exp.template packet<AnyPacket>(index) * AnyPacket(reciprocal(scalar));
        }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return exp.template packetPartial<AnyPacket>(index, count) * AnyPacket(reciprocal(scalar));
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType1, class VectorType2>
    class VectorExpr<ExpressionType::Div, VectorType1, VectorType2>
            : public RValueVector<VectorExpr<ExpressionType::Div, VectorType1, VectorType2>> {
        using This = VectorExpr<ExpressionType::Div, VectorType1, VectorType2>;
        using Base = RValueVector<This>;
        const VectorType1& v1;
        const VectorType2& v2;
    public:
        VectorExpr(const RValueVector<VectorType1>& v1_, const RValueVector<VectorType2>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1.getLength() == v2.getLength());
        }
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return v1.calc(s) / v2.calc(s); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return v1.template packet<AnyPacket>(index) / v2.template packet<AnyPacket>(index);
        }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return v1.template packetPartial<AnyPacket>(index, count) / v2.template packetPartial<AnyPacket>(index, count);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v1.getLength(); }
    };
    //////////////////////////////////////Compare//////////////////////////////////////
    template<class VectorType1, class VectorType2>
    class VectorExpr<ExpressionType::More, VectorType1, VectorType2>
            : public RValueVector<VectorExpr<ExpressionType::More, VectorType1, VectorType2>> {
        using This = VectorExpr<ExpressionType::More, VectorType1, VectorType2>;
        using Base = RValueVector<This>;
        const VectorType1& v1;
        const VectorType2& v2;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpr(const RValueVector<VectorType1>& v1_, const RValueVector<VectorType2>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1.getLength() == v2.getLength());
        }
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(v1.calc(s) > v2.calc(s)); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return v1.template packet<AnyPacket>(index) > v2.template packet<AnyPacket>(index);
        }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return v1.template packetPartial<AnyPacket>(index, count) > v2.template packetPartial<AnyPacket>(index, count);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v1.getLength(); }
    };

    template<class VectorType, class AnyScalar>
    class VectorExpr<ExpressionType::More, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpr<ExpressionType::More, VectorType, ScalarBase<AnyScalar>>> {
        using This = VectorExpr<ExpressionType::More, VectorType, ScalarBase<AnyScalar>>;
        using Base = RValueVector<This>;
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpr(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(exp.calc(s) > scalar); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return exp.template packet<AnyPacket>(index) > AnyPacket(scalar);
        }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return exp.template packetPartial<AnyPacket>(index, count) > AnyPacket(scalar);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType1, class VectorType2>
    class VectorExpr<ExpressionType::MoreEq, VectorType1, VectorType2>
            : public RValueVector<VectorExpr<ExpressionType::MoreEq, VectorType1, VectorType2>> {
        using This = VectorExpr<ExpressionType::MoreEq, VectorType1, VectorType2>;
        using Base = RValueVector<This>;
        const VectorType1& v1;
        const VectorType2& v2;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpr(const RValueVector<VectorType1>& v1_, const RValueVector<VectorType2>& v2_)
                : v1(v1_.getDerived()), v2(v2_.getDerived()) {
            assert(v1.getLength() == v2.getLength());
        }
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(v1.calc(s) >= v2.calc(s)); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return v1.template packet<AnyPacket>(index) >= v2.template packet<AnyPacket>(index);
        }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return v1.template packetPartial<AnyPacket>(index, count) >= v2.template packetPartial<AnyPacket>(index, count);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v1.getLength(); }
    };

    template<class VectorType, class AnyScalar>
    class VectorExpr<ExpressionType::MoreEq, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpr<ExpressionType::MoreEq, VectorType, ScalarBase<AnyScalar>>> {
        using This = VectorExpr<ExpressionType::MoreEq, VectorType, ScalarBase<AnyScalar>>;
        using Base = RValueVector<This>;
        const VectorType& exp;
        const AnyScalar& scalar;
    public:
        using typename Base::ScalarType;
    public:
        VectorExpr(const RValueVector<VectorType>& exp_, const AnyScalar& scalar_)
                : exp(exp_.getDerived()), scalar(scalar_) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(exp.calc(s) >= scalar); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return exp.template packet<AnyPacket>(index) >= AnyPacket(scalar);
        }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return exp.template packetPartial<AnyPacket>(index, count) >= AnyPacket(scalar);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    class VectorExpr<ExpressionType::Reciprocal, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Reciprocal, VectorType>> {
        using This = VectorExpr<ExpressionType::Reciprocal, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& exp;
    public:
        VectorExpr(const RValueVector<VectorType>& exp_) : exp(exp_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t index) const { return reciprocal(exp.calc(index)); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const { return AnyPacket(1) / exp.template packet<AnyPacket>(index); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const { return AnyPacket(1) / exp.template packetPartial<AnyPacket>(index, count); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType>
    class VectorExpr<ExpressionType::Sqrt, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Sqrt, VectorType>> {
        using This = VectorExpr<ExpressionType::Sqrt, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& exp;
    public:
        VectorExpr(const RValueVector<VectorType>& exp_) : exp(exp_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return sqrt(exp.calc(s)); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const { return sqrt(exp.template packet<AnyPacket>(index)); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const { return sqrt(exp.template packetPartial<AnyPacket>(index, count)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType>
    class VectorExpr<ExpressionType::Cbrt, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Cbrt, VectorType>> {
        using This = VectorExpr<ExpressionType::Cbrt, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& exp;
    public:
        VectorExpr(const RValueVector<VectorType>& exp_) : exp(exp_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return cbrt(exp.calc(s)); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            AnyPacket result = exp.template packet<AnyPacket>(index);
            for (size_t i = 0; i < static_cast<size_t>(AnyPacket::size()); ++i)
                result.insert(i, cbrt(result[i]));
            return result;
        }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            AnyPacket result = exp.template packetPartial<AnyPacket>(index, count);
            for (size_t i = 0; i < count; ++i)
                result.insert(i, cbrt(result[i]));
            return result;
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class VectorType>
    class VectorExpr<ExpressionType::Abs, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Abs, VectorType>> {
        using This = VectorExpr<ExpressionType::Abs, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& v;
    public:
        VectorExpr(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return abs(v.calc(s)); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const { return abs(v.template packet<AnyPacket>(index)); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const { return abs(v.template packetPartial<AnyPacket>(index, count)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpr<ExpressionType::Ln, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Ln, VectorType>> {
        using This = VectorExpr<ExpressionType::Ln, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& v;
    public:
        VectorExpr(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return ln(v.calc(s)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpr<ExpressionType::Exp, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Exp, VectorType>> {
        using This = VectorExpr<ExpressionType::Exp, VectorType>;
        using Base = RValueVector<This>;
        const VectorType& v;
    public:
        VectorExpr(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return exp(v.calc(s)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpr<ExpressionType::Pow, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Pow, VectorType>> {
        using This = VectorExpr<ExpressionType::Pow, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
        const ScalarType& s;
    public:
        VectorExpr(const RValueVector<VectorType>& v_, const ScalarBase<ScalarType>& s_)
                : v(v_.getDerived()), s(s_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return pow(v.calc(i), s); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpr<ExpressionType::Sin, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Sin, VectorType>> {
        using This = VectorExpr<ExpressionType::Sin, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
    public:
        VectorExpr(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return sin(v.calc(i)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpr<ExpressionType::Cos, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Cos, VectorType>> {
        using This = VectorExpr<ExpressionType::Cos, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
    public:
        VectorExpr(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return cos(v.calc(i)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpr<ExpressionType::Relu, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Relu, VectorType>> {
        using This = VectorExpr<ExpressionType::Relu, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
    public:
        VectorExpr(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return relu(v.calc(i)); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };

    template<class VectorType>
    class VectorExpr<ExpressionType::Unit, VectorType>
            : public RValueVector<VectorExpr<ExpressionType::Unit, VectorType>> {
        using This = VectorExpr<ExpressionType::Unit, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& v;
    public:
        VectorExpr(const RValueVector<VectorType>& v_) : v(v_.getDerived()) {}
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return v.calc(i).unit(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v.getLength(); }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Div//////////////////////////////////////
    template<class VectorType, class ScalarType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Div, VectorType, ScalarBase<ScalarType>>
    operator/(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) noexcept {
        return VectorExpr<ExpressionType::Div, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline VectorExpr<ExpressionType::Div, VectorType1, VectorType2>
    divide(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpr<ExpressionType::Div, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }
    //////////////////////////////////////Compare//////////////////////////////////////
    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline VectorExpr<ExpressionType::More, VectorType1, VectorType2>
    operator>(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpr<ExpressionType::More, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline VectorExpr<ExpressionType::MoreEq, VectorType1, VectorType2>
    operator>=(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpr<ExpressionType::MoreEq, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    [[nodiscard]] inline VectorExpr<ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>
    operator>=(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) noexcept {
        return VectorExpr<ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline VectorExpr<ExpressionType::More, VectorType2, VectorType1>
    operator<(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpr<ExpressionType::More, VectorType2, VectorType1>(v2.getDerived(), v1.getDerived());
    }

    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline VectorExpr<ExpressionType::MoreEq, VectorType2, VectorType1>
    operator<=(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpr<ExpressionType::MoreEq, VectorType2, VectorType1>(v2.getDerived(), v1.getDerived());
    }

    template<class VectorType, class ScalarType>
    [[nodiscard]] inline VectorExpr<ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>
    operator<=(const ScalarBase<ScalarType>& s, const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::MoreEq, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Reciprocal, VectorType>
    reciprocal(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Reciprocal, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Sqrt, VectorType> sqrt(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Sqrt, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Cbrt, VectorType> cbrt(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Cbrt, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Abs, VectorType> abs(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Abs, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Relu, VectorType> relu(
            const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Relu, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Unit, VectorType> unit(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Unit, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Ln, VectorType> ln(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Ln, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Exp, VectorType> exp(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Exp, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Pow, VectorType> pow(
            const RValueVector<VectorType>& v,
            const ScalarBase<typename VectorType::ScalarType>& s) noexcept {
        return VectorExpr<ExpressionType::Pow, VectorType>(v, s);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Sin, VectorType> sin(
            const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Sin, VectorType>(v);
    }

    template<class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Cos, VectorType> cos(
            const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExpressionType::Cos, VectorType>(v);
    }
}

namespace Physica {
    using namespace Core;

    template<ExpressionType Type, class Expr1, class Expr2>
    class Traits<VectorExpr<Type, Expr1, Expr2>> {
        static_assert(Expr1::SizeAtCompile == Dynamic || Expr2::SizeAtCompile == Dynamic || (Expr1::SizeAtCompile == Expr2::SizeAtCompile)
                         , "[Error]: Vector dimentions do not match");
        using ScalarType1 = typename Expr1::ScalarType;
        using RealType = typename ScalarType1::RealType;
        using BinaryScalarType = typename Core::Internal::BinaryScalarOpReturnType<ScalarType1, typename Expr2::ScalarType>::Type;
        constexpr static bool FastAssign1 = Traits<Expr1>::FastAssign;
        constexpr static bool FastAssign2 = Traits<Expr2>::FastAssign;
        constexpr static bool FastPacket1 = Traits<Expr1>::FastPacket;
        constexpr static bool FastPacket2 = Traits<Expr2>::FastPacket;
        constexpr static bool IsAddOrSub = Type == ExpressionType::Add || Type == ExpressionType::Sub;
    public:
        using ScalarType = typename std::conditional<Type == ExpressionType::Abs, RealType, BinaryScalarType>::type;
        constexpr static size_t SizeAtCompile = Expr1::SizeAtCompile > Expr2::SizeAtCompile ? Expr1::SizeAtCompile : Expr2::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = Expr1::MaxSizeAtCompile > Expr2::MaxSizeAtCompile ? Expr1::MaxSizeAtCompile : Expr2::MaxSizeAtCompile;

        constexpr static bool FastAssign = IsAddOrSub && (FastAssign1 || FastAssign2);
        constexpr static bool FastPacket = FastPacket1 && FastPacket2;
    };

    template<ExpressionType Type, class Expr, class AnyScalar>
    class Traits<VectorExpr<Type, Expr, ScalarBase<AnyScalar>>> {
        static_assert(is_scalar<AnyScalar>::value, "[Error]: This is not a scalar type");
    public:
        using ScalarType = typename Core::Internal::BinaryScalarOpReturnType<typename Expr::ScalarType, AnyScalar>::Type;
        constexpr static size_t SizeAtCompile = Expr::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = Expr::MaxSizeAtCompile;

        constexpr static bool FastAssign = Traits<Expr>::FastAssign;
        constexpr static bool FastPacket = Traits<Expr>::FastPacket;
    };
}

#include "VectorExprImpl/VectorAdd.h"
#include "VectorExprImpl/VectorSub.h"
#include "VectorExprImpl/VectorMul.h"
#include "VectorExprImpl/VectorMinus.h"
#include "VectorExprImpl/VectorSquare.h"
#include "VectorExprImpl/VectorCosh.h"
