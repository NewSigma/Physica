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
    template<Vector V1, Vector V2>
    class VectorExpr<ExprID::More, V1, V2>
            : public BinaryVectorExpr<ExprID::More, V1, V2> {
        using Base = BinaryVectorExpr<ExprID::More, V1, V2>;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t s) const { return T(getLHS().calc(s) > getRHS().calc(s)); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return getLHS().template packet<Pack>(index) > getRHS().template packet<Pack>(index);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return getLHS().template packetPartial<Pack>(index, count) > getRHS().template packetPartial<Pack>(index, count);
        }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector V, Scalar U>
    class VectorExpr<ExprID::More, V, U>
            : public BinaryVectorExpr<ExprID::More, V, U> {
        using Base = BinaryVectorExpr<ExprID::More, V, U>;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t s) const { return T(getLHS().calc(s) > getRHS()); }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return getLHS().template packet<Pack>(index) > Pack(getRHS());
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return getLHS().template packetPartial<Pack>(index, count) > Pack(getRHS());
        }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector V, Scalar U>
    [[nodiscard]] auto operator>(V&& v, U&& x) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprID::More, V&&, U&&>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto operator>(V1&& v1, V2&& v2) noexcept requires(!CUDA<V1> && !CUDA<V2>) {
        return VectorExpr<ExprID::More, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
