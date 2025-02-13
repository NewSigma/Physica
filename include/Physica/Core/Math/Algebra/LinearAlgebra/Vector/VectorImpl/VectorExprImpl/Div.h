/*
 * Copyright 2024-2025 Weibo He.
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

#include "../VectorExpr.h"

namespace Physica {
    template<class T, class U>
    class VectorExpr<ExprType::Div, T, U>
            : public BinaryVectorExpr<ExprType::Div, T, U> {
        static_assert(Scalar<T> || Scalar<U>, "[Error]: Either type should be Scalar");

        using Base = BinaryVectorExpr<ExprType::Div, T, U>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    public:
        VectorExpr(const T& lhs, const U& rhs) : Base(lhs, rhs) {
            if constexpr (Vector<T>)
                assert(!rhs.isZero() && "[Error]: Divide by zero");
        }
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const {
            if constexpr (Vector<T>)
                return Base::getLHS().calc(s) / Base::getRHS();
            else
                return Base::getLHS() / Base::getRHS().calc(s);
        }

        [[nodiscard]] ValueType calc_value(size_t s) const {
            if constexpr (Vector<T>)
                return Base::getLHS().calc_value(s) / Base::getRHS().value();
            else
                return Base::getLHS().value() / Base::getRHS().calc_value(s);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            if constexpr (Vector<T>)
                return Base::getLHS().template packet<AnyPacket>(index) * AnyPacket(reciprocal(Base::getRHS()));
            else
                return AnyPacket(Base::getLHS()) / Base::getRHS().template packet<AnyPacket>(index);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            if constexpr (Vector<T>)
                return Base::getLHS().template packetPartial<AnyPacket>(index, count) * AnyPacket(reciprocal(Base::getRHS()));
            else
                return AnyPacket(Base::getLHS()) / Base::getRHS().template packetPartial<AnyPacket>(index, count);
        }
    };

    template<Vector T1, Vector T2>
    class VectorExpr<ExprType::Div, T1, T2>
            : public BinaryVectorExpr<ExprType::Div, T1, T2> {
        using Base = BinaryVectorExpr<ExprType::Div, T1, T2>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const {
            assert(!Base::getRHS().calc(s).isZero() && "[Error]: Divide by zero");
            return Base::getLHS().calc(s) / Base::getRHS().calc(s);
        }

        [[nodiscard]] ValueType calc_value(size_t s) const {
            assert(!Base::getRHS().calc_value(s).isZero() && "[Error]: Divide by zero");
            return Base::getLHS().calc_value(s) / Base::getRHS().calc_value(s);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return Base::getLHS().template packet<AnyPacket>(index) / Base::getRHS().template packet<AnyPacket>(index);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            const auto pack1 = Base::getLHS().template packetPartial<AnyPacket>(index, count);
            const auto pack2 = Base::getRHS().template packetPartial<AnyPacket>(index, count);
            return (pack1 / pack2).cutoff(count);
        }
    };

    template<Vector T, Scalar U>
    [[nodiscard]] inline auto operator/(const T& v, const U& x) noexcept {
        return VectorExpr<ExprType::Div, T, U>(v, x);
    }

    template<Scalar T, Vector U>
    [[nodiscard]] inline auto operator/(const T& v, const U& x) noexcept {
        return VectorExpr<ExprType::Div, T, U>(v, x);
    }

    template<Vector T1, Vector T2>
    [[nodiscard]] inline auto divide(const T1& v1, const T2& v2) noexcept {
        return VectorExpr<ExprType::Div, T1, T2>(v1, v2);
    }
}
