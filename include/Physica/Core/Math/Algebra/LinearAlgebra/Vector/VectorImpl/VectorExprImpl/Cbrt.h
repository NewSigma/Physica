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
    template<Vector T>
    class VectorExpr<ExprType::Cbrt, T> : public UnitaryVectorExpr<ExprType::Cbrt, T> {
        using This = VectorExpr<ExprType::Cbrt, T>;
        using Base = UnitaryVectorExpr<ExprType::Cbrt, T>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const { return cbrt(Base::getExpr().calc(s)); }

        [[nodiscard]] ValueType calc_value(size_t index) const { return cbrt(Base::getExpr().calc_value(index)); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            Pack result = Base::getExpr().template packet<Pack>(index);
            for (size_t i = 0; i < static_cast<size_t>(Pack::size()); ++i)
                result.insert(i, cbrt(result[i]));
            return result;
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            Pack result = Base::getExpr().template packetPartial<Pack>(index, count);
            for (size_t i = 0; i < count; ++i)
                result.insert(i, cbrt(result[i]));
            return result;
        }
    };

    template<Vector T>
    [[nodiscard]] inline auto cbrt(const T& v) noexcept {
        return VectorExpr<ExprType::Cbrt, T>(v);
    }
}
