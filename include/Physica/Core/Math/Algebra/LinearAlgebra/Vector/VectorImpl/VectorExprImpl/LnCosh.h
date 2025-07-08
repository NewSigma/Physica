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
    template<Vector V>
    class VectorExpr<ExprType::LnCosh, V> : public UnitaryVectorExpr<ExprType::LnCosh, V> {
        using This = VectorExpr<ExprType::LnCosh, V>;
        using Base = UnitaryVectorExpr<ExprType::LnCosh, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const {
            return lncosh(Base::getExpr().calc(index));
        }

        [[nodiscard]] Tv calc_value(size_t index) const {
            return lncosh(Base::getExpr().calc_value(index));
        }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return lncosh(Base::getExpr().template packet<Pack>(index));
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return lncosh(Base::getExpr().template packetPartial<Pack>(index, count));
        }
    };

    template<Vector V>
    [[nodiscard]] auto lncosh(V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::LnCosh, V&&>(std::forward<V>(v));
    }
}
