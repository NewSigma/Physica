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
    class VectorExpr<ExprType::Square, V> : public UnitaryVectorExpr<ExprType::Square, V> {
        using This = VectorExpr<ExprType::Square, V>;
        using Base = UnitaryVectorExpr<ExprType::Square, V>;
    public:
        using Base::isComplex;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Trv;
    public:
        using Base::Base;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& v) const;
        void assign_mkl(Vector auto& v) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_base(Vector auto& v) const;

        [[nodiscard]] CoDiff<T> calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;

        void reverse(const Scalar auto& grad) const noexcept requires(isReverseDiff);
    };

    template<Vector V>
    template<ExecutePolicy P>
    void VectorExpr<ExprType::Square, V>::assign(Vector auto& v) const {
        assign_base<P>(v);
    }

    template<Vector V>
    template<ExecutePolicy P>
    void VectorExpr<ExprType::Square, V>::assign_base(Vector auto& v) const {
        Base::template assign<P>(v);
    }

    template<Vector V>
    auto VectorExpr<ExprType::Square, V>::calc(size_t index) const -> CoDiff<T> {
        return square(Base::getExpr().calc(index));
    }

    template<Vector V>
    auto VectorExpr<ExprType::Square, V>::calc_value(size_t index) const -> Tv {
        return square(Base::getExpr().calc_value(index));
    }

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Square, V>::packet(size_t index) const {
        return square(Base::getExpr().template packet<Pack>(index));
    }

    template<Vector V>
    template<Packet Pack>
    Pack VectorExpr<ExprType::Square, V>::packetPartial(size_t index, size_t count) const {
        return square(Base::getExpr().template packetPartial<Pack>(index, count));
    }

    template<Vector V>
    void VectorExpr<ExprType::Square, V>::reverse(const Scalar auto& grad) const noexcept requires(isReverseDiff) {
        const auto& expr = Base::getExpr();
        expr.reverse(expr.values() * (Trv(2) * grad));
    }

    template<Vector V>
    [[nodiscard]] auto square(V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Square, V&&>(std::forward<V>(v));
    }
}

#ifdef PHYSICA_MKL
    #include "MKL/Square.h"
#endif
