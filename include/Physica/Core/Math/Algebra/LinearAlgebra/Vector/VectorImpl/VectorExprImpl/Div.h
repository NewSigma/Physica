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
    template<Vector T, Scalar U>
    class VectorExpr<ExprType::Div, T, U>
            : public BinaryVectorExpr<ExprType::Div, T, U> {
        using Base = BinaryVectorExpr<ExprType::Div, T, U>;
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

    template<Vector T1, Vector T2>
    class VectorExpr<ExprType::Div, T1, T2>
            : public BinaryVectorExpr<ExprType::Div, T1, T2> {
        using Base = BinaryVectorExpr<ExprType::Div, T1, T2>;
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

    template<Vector T, Scalar U>
    [[nodiscard]] inline auto operator/(const T& v, const U& x) noexcept {
        return VectorExpr<ExprType::Div, T, U>(v, x);
    }

    template<Vector T1, Vector T2>
    [[nodiscard]] inline auto divide(const T1& v1, const T2& v2) noexcept {
        return VectorExpr<ExprType::Div, T1, T2>(v1, v2);
    }
}
