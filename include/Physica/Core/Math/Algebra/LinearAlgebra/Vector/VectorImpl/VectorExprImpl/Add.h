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
    template<Vector T, Scalar U>
    class VectorExpr<ExprType::Add, T, U>
            : public BinaryVectorExpr<ExprType::Add, T, U> {
        using Base = BinaryVectorExpr<ExprType::Add, T, U>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const {
            return Base::getLHS().calc(s) + Base::getRHS();
        }

        [[nodiscard]] ValueType calc_value(size_t s) const {
            return Base::getLHS().calc_value(s) + Base::getRHS().value();
        }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return Base::getLHS().template packet<Pack>(index) + Pack(Base::getRHS());
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<Pack>(index, count) + Pack(Base::getRHS(), count);
        }

        template<Vector V>
        void reverse(const V& grad_) const noexcept requires(isReverseDiff) {
            const auto& grad = grad_.values();
            if constexpr (ReverseDiff<T>)
                Base::getLHS().reverse(grad);
            if constexpr (ReverseDiff<U>)
                Base::getRHS().reverse(grad.sum());
        }
    };

    template<Vector T1, Vector T2>
    class VectorExpr<ExprType::Add, T1, T2>
            : public BinaryVectorExpr<ExprType::Add, T1, T2> {
        using Base = BinaryVectorExpr<ExprType::Add, T1, T2>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        template<Vector V, class Executor = SeqExecutor>
        inline void assign(V& v) const;

        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const {
            return getLHS().calc(s) + getRHS().calc(s);
        }

        [[nodiscard]] ValueType calc_value(size_t s) const {
            return Base::getLHS().calc_value(s) + Base::getRHS().calc_value(s);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return getLHS().template packet<Pack>(index) + getRHS().template packet<Pack>(index);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return getLHS().template packetPartial<Pack>(index, count) + getRHS().template packetPartial<Pack>(index, count);
        }

        template<class U>
        void reverse(const U& grad_) const noexcept requires(isReverseDiff);
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector T1, Vector T2>
    template<Vector V, class Executor>
    inline void VectorExpr<ExprType::Add, T1, T2>::assign(V& v) const {
        constexpr bool FastAssign1 = Traits<T1>::FastAssign;
        constexpr bool FastAssign2 = Traits<T2>::FastAssign;
        if constexpr (FastAssign1) {
            getLHS().template assign<V, Executor>(v);
            v += getRHS();
        }
        else if constexpr (FastAssign2) {
                getRHS().template assign<V, Executor>(v);
                v += getLHS();
        }
        else
            Base::template assign<V, Executor>(v);
    }

    template<Vector T1, Vector T2>
    template<class U>
    void VectorExpr<ExprType::Add, T1, T2>::reverse(const U& grad_) const noexcept requires(isReverseDiff) {
        if constexpr (Scalar<U>) {
            if constexpr (ReverseDiff<T1>)
                Base::getLHS().reverse(grad_);
            if constexpr (ReverseDiff<T2>)
                Base::getRHS().reverse(grad_);
        }
        else {
            static_assert(Vector<U>);
            const auto& grad = grad_.values();
            if constexpr (ReverseDiff<T1>)
                Base::getLHS().reverse(grad);
            if constexpr (ReverseDiff<T2>)
                Base::getRHS().reverse(grad);
        }
    }

    template<Vector T, Scalar U>
    [[nodiscard]] inline auto operator+(const T& v, const U& x) noexcept {
        return VectorExpr<ExprType::Add, T, U>(v, x);
    }

    template<Scalar U, Vector T>
    [[nodiscard]] inline auto operator+(const U& x, const T& v) noexcept {
        return v + x;
    }

    template<Vector T1, Vector T2>
    [[nodiscard]] inline auto operator+(const T1& v1, const T2& v2) noexcept {
        return VectorExpr<ExprType::Add, T1, T2>(v1, v2);
    }
}
