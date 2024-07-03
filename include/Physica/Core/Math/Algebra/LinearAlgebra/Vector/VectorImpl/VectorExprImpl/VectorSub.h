/*
 * Copyright 2024 WeiBo He.
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
    class VectorExpr<ExpressionType::Sub, VectorType1, VectorType2>
            : public RValueVector<VectorExpr<ExpressionType::Sub, VectorType1, VectorType2>> {
        using This = VectorExpr<ExpressionType::Sub, VectorType1, VectorType2>;
    public:
        using Base = RValueVector<This>;
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
        inline void assignTo(LValueVector<OtherDerived>& v_) const;

        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return ScalarType(expr1.calc(s)) - ScalarType(expr2.calc(s)); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return expr1.template packet<PacketType>(index) - expr2.template packet<PacketType>(index);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const {
            return expr1.template packetPartial<PacketType>(index, count) - expr2.template packetPartial<PacketType>(index, count);
        }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const { return expr1.getLength(); }
    };

    template<class VectorType1, class VectorType2>
    template<class OtherDerived, class Executor>
    inline void VectorExpr<ExpressionType::Sub, VectorType1, VectorType2>::assignTo(LValueVector<OtherDerived>& v_) const {
        constexpr bool FastAssign1 = Traits<VectorType1>::FastAssign;
        constexpr bool FastAssign2 = Traits<VectorType2>::FastAssign;
        auto& v = v_.getDerived();
        if constexpr (FastAssign1) {
            if constexpr (FastAssign2) {
                expr1.template assignTo<OtherDerived, Executor>(v);
                OtherDerived buffer;
                buffer.template operator=<VectorType2, Executor>(expr2);
                v -= buffer;
            }
            else {
                expr1.template assignTo<OtherDerived, Executor>(v);
                v -= expr2;
            }
        }
        else {
            if constexpr (FastAssign2) {
                (-expr2).template assignTo<OtherDerived, Executor>(v);
                v += expr1;
            }
            else
                Base::template assignTo<OtherDerived, Executor>(v);
        }
    }

    template<class VectorType, class AnyScalar>
    class VectorExpr<ExpressionType::Sub, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpr<ExpressionType::Sub, VectorType, ScalarBase<AnyScalar>>> {
        using This = VectorExpr<ExpressionType::Sub, VectorType, ScalarBase<AnyScalar>>;
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
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return exp.calc(s) - scalar; }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return exp.template packet<PacketType>(index) - PacketType(scalar);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const {
            return exp.template packetPartial<PacketType>(index, count) - PacketType(scalar);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return exp.getLength(); }
    };

    template<class Derived, class OtherDerived>
    [[nodiscard]] inline VectorExpr<ExpressionType::Sub, Derived, OtherDerived>
    operator-(const RValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) noexcept {
        return VectorExpr<ExpressionType::Sub, Derived, OtherDerived>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    [[nodiscard]] inline VectorExpr<ExpressionType::Sub, VectorType, ScalarBase<ScalarType>>
    operator-(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) noexcept {
        return VectorExpr<ExpressionType::Sub, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }
}
