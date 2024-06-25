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
    template<class VectorType, class AnyScalar>
    class VectorExpression<ExpressionType::Mul, VectorType, ScalarBase<AnyScalar>>
            : public RValueVector<VectorExpression<ExpressionType::Mul, VectorType, ScalarBase<AnyScalar>>> {
        using This = VectorExpression<ExpressionType::Mul, VectorType, ScalarBase<AnyScalar>>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const VectorType& expr;
        const AnyScalar& s;
    public:
        VectorExpression(const RValueVector<VectorType>& expr_, const AnyScalar& s_)
                : expr(expr_.getDerived()), s(s_) {}
        VectorExpression(const This&) = delete;
        VectorExpression(This&&) noexcept = delete;
        ~VectorExpression() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherDerived, class Executor = SequentialExecutor>
        inline void assignTo(LValueVector<OtherDerived>& v) const;

        [[nodiscard]] ScalarType calc(size_t index) const { return ScalarType(expr.calc(index)) * ScalarType(s); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return expr.template packet<PacketType>(index) * PacketType(s);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const {
            return expr.template packetPartial<PacketType>(index, count) * PacketType(s);
        }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const { return expr.getLength(); }
    };

    template<class VectorType, class AnyScalar>
    template<class OtherDerived, class Executor>
    inline void VectorExpression<ExpressionType::Mul, VectorType, ScalarBase<AnyScalar>>::assignTo(LValueVector<OtherDerived>& v) const {
        constexpr bool FastAssign = Internal::Traits<VectorType>::FastAssign;
        if constexpr (FastAssign) {
            expr.assignTo(v.getDerived());
            v.getDerived() *= s;
        }
        else
            Base::assignTo(v);
    }

    template<class VectorType1, class VectorType2>
    class VectorExpression<ExpressionType::Mul, VectorType1, VectorType2>
            : public RValueVector<VectorExpression<ExpressionType::Mul, VectorType1, VectorType2>> {
        using This = VectorExpression<ExpressionType::Mul, VectorType1, VectorType2>;
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
        [[nodiscard]] typename Base::ScalarType calc(size_t index) const { return v1.calc(index) * v2.calc(index); }
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const {
            return v1.template packet<PacketType>(index) * v2.template packet<PacketType>(index);
        }
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index, size_t count) const {
            return v1.template packetPartial<PacketType>(index, count) * v2.template packetPartial<PacketType>(index, count);
        }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const { return v1.getLength(); }
    };

    template<class VectorType, class ScalarType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Mul, VectorType, ScalarBase<ScalarType>>
    operator*(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) noexcept {
        return VectorExpression<ExpressionType::Mul, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class ScalarType, class VectorType>
    [[nodiscard]] inline VectorExpression<ExpressionType::Mul, VectorType, ScalarBase<ScalarType>>
    operator*(const ScalarBase<ScalarType>& s, const RValueVector<VectorType>& v) noexcept {
        return v * s;
    }
    
    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline VectorExpression<ExpressionType::Mul, VectorType1, VectorType2> hadamard(
            const RValueVector<VectorType1>& v1,
            const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpression<ExpressionType::Mul, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }
}
