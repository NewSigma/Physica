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
    class VectorExpr<ExpressionType::Minus, VectorType>
            : public UnitaryVectorExpr<ExpressionType::Minus, VectorType> {
        using Base = UnitaryVectorExpr<ExpressionType::Minus, VectorType>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] typename Base::ScalarType calc(size_t s) const { return -Base::getExpr().calc(s); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const { return -Base::getExpr().template packet<AnyPacket>(index); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const { return -Base::getExpr().template packetPartial<AnyPacket>(index, count); }
    };

    template<class Derived>
    [[nodiscard]] inline auto operator-(const RValueVector<Derived>& v) noexcept {
        return VectorExpr<ExpressionType::Minus, Derived>(v.getDerived());
    }
}
