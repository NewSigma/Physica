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
    class VectorExpr<ExprType::Minus, T>
            : public UnitaryVectorExpr<ExprType::Minus, T> {
        using Base = UnitaryVectorExpr<ExprType::Minus, T>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const { return -Base::getExpr().calc(s); }

        [[nodiscard]] ValueType calc_value(size_t index) const { return -Base::getExpr().calc_value(index); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const { return -Base::getExpr().template packet<Pack>(index); }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const { return -Base::getExpr().template packetPartial<Pack>(index, count); }

        template<class U>
        void reverse(const U& grad_) const noexcept requires(isReverseDiff) {
            Base::getExpr().reverse(-grad_);
        }
    };

    template<Vector T>
    [[nodiscard]] inline auto operator-(const T& v) noexcept {
        return VectorExpr<ExprType::Minus, T>(v);
    }
}
