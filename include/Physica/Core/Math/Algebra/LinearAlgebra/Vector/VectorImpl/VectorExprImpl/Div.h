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
    template<class U, class V>
    class VectorExpr<ExprType::Div, U, V>
            : public BinaryVectorExpr<ExprType::Div, U, V> {
        static_assert(Scalar<U> || Scalar<V>, "[Error]: Either type should be Scalar");

        using Base = BinaryVectorExpr<ExprType::Div, U, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        VectorExpr(U lhs, V rhs);
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t s) const;
        [[nodiscard]] Tv calc_value(size_t s) const;

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;
    };

    template<class U, class V>
    VectorExpr<ExprType::Div, U, V>::VectorExpr(U lhs, V rhs) : Base(std::forward<U>(lhs), std::forward<V>(rhs)) {
        if constexpr (Vector<U>)
            assert(!Base::getRHS().isZero() && "[Error]: Divide by zero");
    }

    template<class U, class V>
    auto VectorExpr<ExprType::Div, U, V>::calc(size_t s) const -> CoDiff<T> {
        if constexpr (Vector<U>)
            return Base::getLHS().calc(s) / Base::getRHS();
        else
            return Base::getLHS() / Base::getRHS().calc(s);
    }

    template<class U, class V>
    auto VectorExpr<ExprType::Div, U, V>::calc_value(size_t s) const -> Tv {
        if constexpr (Vector<U>)
            return Base::getLHS().calc_value(s) / Base::getRHS().value();
        else
            return Base::getLHS().value() / Base::getRHS().calc_value(s);
    }

    template<class U, class V>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Div, U, V>::packet(size_t index) const {
        if constexpr (Vector<U>)
            return Base::getLHS().template packet<Pack>(index) * Pack(reciprocal(Base::getRHS()));
        else {
            auto div = Base::getRHS().template packet<Pack>(index);
            assert(!div.isZero().horizontal_or() && "[Error]: Divide by zero");
            return Pack(Base::getLHS()) / div;
        }
    }

    template<class U, class V>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Div, U, V>::packetPartial(size_t index, size_t count) const {
        if constexpr (Vector<U>)
            return Base::getLHS().template packetPartial<Pack>(index, count) * Pack(reciprocal(Base::getRHS()));
        else {
            auto div = Base::getRHS().template packetPartial<Pack>(index, count);
            for (size_t i = 0; i < count; ++i)
                assert(!div[i].isZero() && "[Error]: Divide by zero");
            return (Pack(Base::getLHS()) / div).cutoff(count);
        }
    }

    template<Vector V1, Vector V2>
    class VectorExpr<ExprType::Div, V1, V2>
            : public BinaryVectorExpr<ExprType::Div, V1, V2> {
        using Base = BinaryVectorExpr<ExprType::Div, V1, V2>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t s) const;
        [[nodiscard]] Tv calc_value(size_t s) const;

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;
    };

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprType::Div, V1, V2>::calc(size_t s) const -> CoDiff<T> {
        assert(!Base::getRHS().calc(s).isZero() && "[Error]: Divide by zero");
        return Base::getLHS().calc(s) / Base::getRHS().calc(s);
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprType::Div, V1, V2>::calc_value(size_t s) const -> Tv {
        assert(!Base::getRHS().calc_value(s).isZero() && "[Error]: Divide by zero");
        return Base::getLHS().calc_value(s) / Base::getRHS().calc_value(s);
    }

    template<Vector V1, Vector V2>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Div, V1, V2>::packet(size_t index) const {
        auto div = Base::getRHS().template packet<Pack>(index);
        assert(!div.isZero().horizontal_or() && "[Error]: Divide by zero");
        return Base::getLHS().template packet<Pack>(index) / div;
    }

    template<Vector V1, Vector V2>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Div, V1, V2>::packetPartial(size_t index, size_t count) const {
        const auto pack1 = Base::getLHS().template packetPartial<Pack>(index, count);
        const auto pack2 = Base::getRHS().template packetPartial<Pack>(index, count);
        for (size_t i = 0; i < count; ++i)
            assert(!pack2[i].isZero() && "[Error]: Divide by zero");
        return (pack1 / pack2).cutoff(count);
    }

    template<Vector V, Scalar U>
    [[nodiscard]] auto operator/(V&& v, U&& x) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Div, V&&, U&&>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard]] auto operator/(U&& x, V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Div, U&&, V&&>(std::forward<U>(x), std::forward<V>(v));
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto divide(V1&& v1, V2&& v2) noexcept requires(!CUDA<V1> && !CUDA<V2>) {
        return VectorExpr<ExprType::Div, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
