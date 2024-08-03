/*
 * Copyright 2024 Weibo He.
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

namespace Physica::Core {
    template<class VectorType1, class VectorType2>
    class VectorExpr<ExpressionType::Add, VectorType1, VectorType2>
            : public RValueVector<VectorExpr<ExpressionType::Add, VectorType1, VectorType2>> {
        using This = VectorExpr<ExpressionType::Add, VectorType1, VectorType2>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType1& expr1;
        const VectorType2& expr2;
    public:
        VectorExpr(const RValueVector<VectorType1>& expr1_, const RValueVector<VectorType2>& expr2_)
                : expr1(expr1_.getDerived()), expr2(expr2_.getDerived()) {
            assert(expr1.getLength() == expr2.getLength());
        }
        VectorExpr(const This&) = delete;
        VectorExpr(This&&) noexcept = delete;
        ~VectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherDerived, class Executor = SequentialExecutor>
        inline void assignTo(LValueVector<OtherDerived>& v) const;

        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(expr1.calc(s)) + ScalarType(expr2.calc(s)); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return expr1.template packet<AnyPacket>(index) + expr2.template packet<AnyPacket>(index);
        }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return expr1.template packetPartial<AnyPacket>(index, count) + expr2.template packetPartial<AnyPacket>(index, count);
        }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const { return expr1.getLength(); }
    };

    template<class VectorType1, class VectorType2>
    template<class OtherDerived, class Executor>
    inline void VectorExpr<ExpressionType::Add, VectorType1, VectorType2>::assignTo(LValueVector<OtherDerived>& v) const {
        constexpr bool FastAssign1 = Traits<VectorType1>::FastAssign;
        constexpr bool FastAssign2 = Traits<VectorType2>::FastAssign;
        if constexpr (FastAssign1) {
            if constexpr (FastAssign2) {
                expr1.template assignTo<OtherDerived, Executor>(v);
                OtherDerived buffer;
                buffer.template operator=<VectorType2, Executor>(expr2);
                v += buffer;
            }
            else {
                expr1.template assignTo<OtherDerived, Executor>(v);
                v += expr2;
            }
        }
        else {
            if constexpr (FastAssign2) {
                expr2.template assignTo<OtherDerived, Executor>(v);
                v += expr1;
            }
            else
                Base::template assignTo<OtherDerived, Executor>(v);
        }
    }

    template<class VectorType, class AnyScalar>
    class VectorExpr<ExpressionType::Add, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpr<ExpressionType::Add, VectorType, ScalarBase<AnyScalar>>> {
        using This = VectorExpr<ExpressionType::Add, VectorType, ScalarBase<AnyScalar>>;
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
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(exp.calc(s)) + ScalarType(scalar); }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return exp.template packet<AnyPacket>(index) + AnyPacket(scalar);
        }
        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return exp.template packetPartial<AnyPacket>(index, count) + AnyPacket(scalar);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class Derived, class OtherDerived>
    [[nodiscard]] inline VectorExpr<ExpressionType::Add, Derived, OtherDerived>
    operator+(const RValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) noexcept {
        return VectorExpr<ExpressionType::Add, Derived, OtherDerived>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Add, VectorType, ScalarBase<ScalarType>>
    operator+(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) noexcept {
        return VectorExpr<ExpressionType::Add, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class ScalarType, class VectorType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Add, VectorType, ScalarBase<ScalarType>>
    operator+(const ScalarBase<ScalarType>& s, const RValueVector<VectorType>& v) noexcept {
        return v + s;
    }
}
