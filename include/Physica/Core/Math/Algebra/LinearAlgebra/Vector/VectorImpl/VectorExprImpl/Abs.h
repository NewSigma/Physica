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
    template<class VectorType>
    class VectorExpr<ExprType::Abs, VectorType> : public UnitaryVectorExpr<ExprType::Abs, VectorType> {
        using This = VectorExpr<ExprType::Abs, VectorType>;
        using Base = UnitaryVectorExpr<ExprType::Abs, VectorType>;
    public:
        using typename Base::ScalarType;
    private:
        constexpr static bool isComplexV = VectorType::isComplex;
        constexpr static bool isReverseDiff = Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const { return abs(Base::getExpr().calc(index)); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            if constexpr (isComplexV)
                return sqrt(toSquaredNormVector(Base::getExpr()).template packet<AnyPacket>(index));
            else
                return abs(Base::getExpr().template packet<AnyPacket>(index));
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            if constexpr (isComplexV)
                return sqrt(toSquaredNormVector(Base::getExpr()).template packetPartial<AnyPacket>(index, count));
            else
                return abs(Base::getExpr().template packetPartial<AnyPacket>(index, count));
        }
    };

    template<class VectorType>
    [[nodiscard]] inline auto abs(const RValueVector<VectorType>& v) noexcept {
        return VectorExpr<ExprType::Abs, VectorType>(v);
    }
}
