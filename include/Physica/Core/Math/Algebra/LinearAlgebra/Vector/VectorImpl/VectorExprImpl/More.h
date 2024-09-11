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
    class VectorExpr<ExprType::More, VectorType1, VectorType2>
            : public BinaryVectorExpr<ExprType::More, VectorType1, VectorType2> {
        using Base = BinaryVectorExpr<ExprType::More, VectorType1, VectorType2>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(getLHS().calc(s) > getRHS().calc(s)); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return getLHS().template packet<AnyPacket>(index) > getRHS().template packet<AnyPacket>(index);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return getLHS().template packetPartial<AnyPacket>(index, count) > getRHS().template packetPartial<AnyPacket>(index, count);
        }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<class VectorType, class AnyScalar>
    class VectorExpr<ExprType::More, VectorType, ScalarBase<AnyScalar>>
            : public BinaryVectorExpr<ExprType::More, VectorType, ScalarBase<AnyScalar>> {
        using Base = BinaryVectorExpr<ExprType::More, VectorType, ScalarBase<AnyScalar>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(getLHS().calc(s) > getRHS()); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return getLHS().template packet<AnyPacket>(index) > AnyPacket(getRHS());
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return getLHS().template packetPartial<AnyPacket>(index, count) > AnyPacket(getRHS());
        }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline auto operator>(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpr<ExprType::More, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    [[nodiscard]] inline auto operator>(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) noexcept {
        return VectorExpr<ExprType::More, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }
}
