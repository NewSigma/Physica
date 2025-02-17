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

namespace Physica {
    template<Vector T1, Vector T2>
    class VectorExpr<ExprType::MoreEq, T1, T2>
            : public BinaryVectorExpr<ExprType::MoreEq, T1, T2> {
        using Base = BinaryVectorExpr<ExprType::MoreEq, T1, T2>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(getLHS().calc(s) >= getRHS().calc(s)); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return getLHS().template packet<Pack>(index) >= getRHS().template packet<Pack>(index);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return getLHS().template packetPartial<Pack>(index, count) >= getRHS().template packetPartial<Pack>(index, count);
        }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector T, Scalar U>
    class VectorExpr<ExprType::MoreEq, T, U>
            : public BinaryVectorExpr<ExprType::MoreEq, T, U> {
        using Base = BinaryVectorExpr<ExprType::MoreEq, T, U>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t s) const { return ScalarType(getLHS().calc(s) >= getRHS()); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return getLHS().template packet<Pack>(index) >= Pack(getRHS());
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return getLHS().template packetPartial<Pack>(index, count) >= Pack(getRHS());
        }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector T1, Vector T2>
    [[nodiscard]] inline auto operator>=(const T1& v1, const T2& v2) noexcept {
        return VectorExpr<ExprType::MoreEq, T1, T2>(v1, v2);
    }

    template<Vector T, Scalar U>
    [[nodiscard]] inline auto operator>=(const T& v, const U& x) noexcept {
        return VectorExpr<ExprType::MoreEq, T, U>(v, x);
    }
}
