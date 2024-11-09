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
    template<class VectorType, class AnyScalar>
    class VectorExpr<ExprType::Mul, VectorType, ScalarBase<AnyScalar>>
            : public BinaryVectorExpr<ExprType::Mul, VectorType, ScalarBase<AnyScalar>> {
        using Base = BinaryVectorExpr<ExprType::Mul, VectorType, ScalarBase<AnyScalar>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        template<class OtherDerived, class Executor = SequentialExecutor>
        inline void assignTo(LValueVector<OtherDerived>& v) const;

        [[nodiscard]] ScalarType calc(size_t index) const { return ScalarType(Base::getLHS().calc(index)) * ScalarType(Base::getRHS()); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return Base::getLHS().template packet<AnyPacket>(index) * Base::getRHS();
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<AnyPacket>(index, count) * Base::getRHS();
        }

        [[nodiscard]] ScalarType sum() const { return Base::getLHS().sum() * Base::getRHS(); }
    };

    template<class VectorType, class AnyScalar>
    template<class OtherDerived, class Executor>
    inline void VectorExpr<ExprType::Mul, VectorType, ScalarBase<AnyScalar>>::assignTo(LValueVector<OtherDerived>& v) const {
        constexpr bool FastAssign = Traits<VectorType>::FastAssign;
        if constexpr (FastAssign) {
            Base::getLHS().template assignTo<OtherDerived, Executor>(v.getDerived());
            v.getDerived() *= Base::getRHS();
        }
        else
            Base::assignTo(v);
    }

    template<class VectorType1, class VectorType2>
    class VectorExpr<ExprType::Mul, VectorType1, VectorType2>
            : public BinaryVectorExpr<ExprType::Mul, VectorType1, VectorType2> {
        using Base = BinaryVectorExpr<ExprType::Mul, VectorType1, VectorType2>;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t index) const {
            return Base::getLHS().calc(index) * Base::getRHS().calc(index);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return Base::getLHS().template packet<AnyPacket>(index) * Base::getRHS().template packet<AnyPacket>(index);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<AnyPacket>(index, count) * Base::getRHS().template packetPartial<AnyPacket>(index, count);
        }
    };

    template<class VectorType, class ScalarType>
    [[nodiscard]] inline auto operator*(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) noexcept {
        return VectorExpr<ExprType::Mul, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class ScalarType, class VectorType>
    [[nodiscard]] inline auto operator*(const ScalarBase<ScalarType>& s, const RValueVector<VectorType>& v) noexcept {
        return v * s;
    }
    
    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline auto hadamard(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpr<ExprType::Mul, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }
}
