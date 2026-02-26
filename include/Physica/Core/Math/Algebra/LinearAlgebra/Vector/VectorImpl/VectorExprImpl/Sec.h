/*
 * Copyright 2024-2026 Weibo He.
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
    template<Vector V>
    class VectorExpr<ExprID::Sec, V> : public UnitaryVectorExpr<ExprID::Sec, V> {
        using This = VectorExpr<ExprID::Sec, V>;
        using Base = UnitaryVectorExpr<ExprID::Sec, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const {
            return sec(Base::getExpr().calc(index));
        }

        [[nodiscard]] Tv calc_value(size_t index) const {
            return sec(Base::getExpr().calc_value(index));
        }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const noexcept { return sec(Base::getExpr().template packet<Pack>(index)); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index, size_t count) const noexcept {
            return sec(Base::getExpr().template packet<Pack>(index, count)).cutoff(count);
        }
    };

    template<Vector V>
    [[nodiscard]] auto sec(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Sec, V&&>(std::forward<V>(v));
    }
}
