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

#include "../VectorExpr.h"

namespace Physica::Core {
    template<Vector T, Scalar U>
    class VectorExpr<ExprType::Add, T, U>
            : public BinaryVectorExpr<ExprType::Add, T, U> {
        using Base = BinaryVectorExpr<ExprType::Add, T, U>;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const { return ScalarType(Base::getLHS().calc(s)) + ScalarType(Base::getRHS()); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return Base::getLHS().template packet<AnyPacket>(index) + AnyPacket(Base::getRHS());
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<AnyPacket>(index, count) + AnyPacket(Base::getRHS(), count);
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
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        template<LVector V, class Executor = SequentialExecutor>
        inline void assignTo(V& v) const;

        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const { return ScalarType(getLHS().calc(s)) + ScalarType(getRHS().calc(s)); }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return getLHS().template packet<AnyPacket>(index) + getRHS().template packet<AnyPacket>(index);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return getLHS().template packetPartial<AnyPacket>(index, count) + getRHS().template packetPartial<AnyPacket>(index, count);
        }

        template<Vector V>
        void reverse(const V& grad_) const noexcept requires(isReverseDiff) {
            const auto& grad = grad_.values();
            if constexpr (ReverseDiff<T1>)
                Base::getLHS().reverse(grad);
            if constexpr (ReverseDiff<T2>)
                Base::getRHS().reverse(grad);
        }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector T1, Vector T2>
    template<LVector V, class Executor>
    inline void VectorExpr<ExprType::Add, T1, T2>::assignTo(V& v) const {
        constexpr bool FastAssign1 = Traits<T1>::FastAssign;
        constexpr bool FastAssign2 = Traits<T2>::FastAssign;
        if constexpr (FastAssign1) {
            if constexpr (FastAssign2) {
                getLHS().template assignTo<V, Executor>(v);
                V buffer;
                buffer.template operator=<T2, Executor>(getRHS());
                v += buffer;
            }
            else {
                getLHS().template assignTo<V, Executor>(v);
                v += getRHS();
            }
        }
        else {
            if constexpr (FastAssign2) {
                getRHS().template assignTo<V, Executor>(v);
                v += getLHS();
            }
            else
                Base::template assignTo<V, Executor>(v);
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
