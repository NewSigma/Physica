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
    class VectorExpr<ExprType::Tan, T> : public UnitaryVectorExpr<ExprType::Tan, T> {
        using This = VectorExpr<ExprType::Tan, T>;
        using Base = UnitaryVectorExpr<ExprType::Tan, T>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t index) const { return tan(Base::getExpr().calc(index)); }

        [[nodiscard]] ValueType calc_value(size_t index) const { return tan(Base::getExpr().calc_value(index)); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const { return tan(Base::getExpr().template packet<Pack>(index)); }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const { return tan(Base::getExpr().template packetPartial<Pack>(index, count)); }
    };

    template<Vector T>
    [[nodiscard]] inline auto tan(T&& v) noexcept requires(!CUDA<T>) {
        return VectorExpr<ExprType::Tan, T&&>(std::forward<T>(v));
    }
}
