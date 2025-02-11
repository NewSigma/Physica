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
    template<Vector T>
    class VectorExpr<ExprType::Tanh, T> : public UnitaryVectorExpr<ExprType::Tanh, T> {
        using This = VectorExpr<ExprType::Tanh, T>;
        using Base = UnitaryVectorExpr<ExprType::Tanh, T>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t index) const { return tanh(Base::getExpr().calc(index)); }

        [[nodiscard]] ValueType calc_value(size_t index) const { return tanh(Base::getExpr().calc_value(index)); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return tanh(Base::getExpr().template packet<AnyPacket>(index));
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return tanh(Base::getExpr().template packetPartial<AnyPacket>(index, count));
        }

        template<Vector V>
        void reverse(const V& grad_) const noexcept requires(isReverseDiff) {
            const auto& expr = Base::getExpr();
            expr.reverse(hadamard(square(sech(expr.values())), grad_));
        }
    };

    template<Vector T>
    [[nodiscard]] inline auto tanh(const T& v) noexcept {
        return VectorExpr<ExprType::Tanh, T>(v);
    }
}
