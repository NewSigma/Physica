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
    class VectorExpr<ExprType::Div, VectorType, ScalarBase<AnyScalar>>
            : public BinaryVectorExpr<ExprType::Div, VectorType, ScalarBase<AnyScalar>> {
        using Base = BinaryVectorExpr<ExprType::Div, VectorType, ScalarBase<AnyScalar>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const {
            assert(!Base::getRHS().isZero() && "[Error]: Divide by zero");
            return Base::getLHS().calc(s) / Base::getRHS();
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return Base::getLHS().template packet<AnyPacket>(index) * AnyPacket(reciprocal(Base::getRHS()));
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<AnyPacket>(index, count) * AnyPacket(reciprocal(Base::getRHS()));
        }
    };

    template<class VectorType1, class VectorType2>
    class VectorExpr<ExprType::Div, VectorType1, VectorType2>
            : public BinaryVectorExpr<ExprType::Div, VectorType1, VectorType2> {
        using Base = BinaryVectorExpr<ExprType::Div, VectorType1, VectorType2>;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return Base::getLHS().calc(s) / Base::getRHS().calc(s); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return Base::getLHS().template packet<AnyPacket>(index) / Base::getRHS().template packet<AnyPacket>(index);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<AnyPacket>(index, count) / Base::getRHS().template packetPartial<AnyPacket>(index, count);
        }
    };

    template<class VectorType, class ScalarType>
    [[nodiscard]] inline auto operator/(const RValueVector<VectorType>& v, const ScalarBase<ScalarType>& s) noexcept {
        return VectorExpr<ExprType::Div, VectorType, ScalarBase<ScalarType>>(v.getDerived(), s.getDerived());
    }

    template<class VectorType1, class VectorType2>
    [[nodiscard]] inline auto divide(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) noexcept {
        return VectorExpr<ExprType::Div, VectorType1, VectorType2>(v1.getDerived(), v2.getDerived());
    }
}
