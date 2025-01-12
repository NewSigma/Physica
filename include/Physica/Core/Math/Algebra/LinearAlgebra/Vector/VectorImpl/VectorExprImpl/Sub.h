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

namespace Physica::Core {
    template<class T, class U>
    class VectorExpr<ExprType::Sub, T, U>
            : public BinaryVectorExpr<ExprType::Sub, T, U> {
        static_assert(Scalar<T> || Scalar<U>, "[Error]: Either type should be Scalar");

        using Base = BinaryVectorExpr<ExprType::Sub, T, U>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const {
            if constexpr (Vector<T>)
                return Base::getLHS().calc(s) - Base::getRHS();
            else
                return Base::getLHS() - Base::getRHS().calc(s);
        }

        [[nodiscard]] ValueType calc_value(size_t s) const {
            if constexpr (Vector<T>)
                return Base::getLHS().calc_value(s) - Base::getRHS().value();
            else
                return Base::getLHS().value() - Base::getRHS().calc_value(s);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            if constexpr (Vector<T>)
                return Base::getLHS().template packet<AnyPacket>(index) - AnyPacket(Base::getRHS());
            else
                return AnyPacket(Base::getLHS()) - Base::getRHS().template packet<AnyPacket>(index);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            if constexpr (Vector<T>)
                return Base::getLHS().template packetPartial<AnyPacket>(index, count) - AnyPacket(Base::getRHS(), count);
            else
                return AnyPacket(Base::getLHS(), count) - Base::getRHS().template packetPartial<AnyPacket>(index, count);
        }

        template<Vector V>
        void reverse(const V& grad_) const noexcept requires(isReverseDiff) {
            const auto& grad = grad_.values();
            if constexpr (Vector<T>) {
                if constexpr (ReverseDiff<T>)
                    Base::getLHS().reverse(grad);
                if constexpr (ReverseDiff<U>)
                    Base::getRHS().reverse(-grad.sum());
            }
            else {
                if constexpr (ReverseDiff<T>)
                    Base::getLHS().reverse(grad.sum());
                if constexpr (ReverseDiff<U>)
                    Base::getRHS().reverse(-grad);
            }
        }
    };

    template<Vector T1, Vector T2>
    class VectorExpr<ExprType::Sub, T1, T2>
            : public BinaryVectorExpr<ExprType::Sub, T1, T2> {
        using Base = BinaryVectorExpr<ExprType::Sub, T1, T2>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        template<Vector V, class Executor = SequentialExecutor>
        inline void assign(V& v) const;

        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const {
            return getLHS().calc(s) - getRHS().calc(s);
        }

        [[nodiscard]] ValueType calc_value(size_t s) const {
            return getLHS().calc_value(s) - getRHS().calc_value(s);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packet(size_t index) const {
            return getLHS().template packet<AnyPacket>(index) - getRHS().template packet<AnyPacket>(index);
        }

        template<class AnyPacket>
        [[nodiscard]] AnyPacket packetPartial(size_t index, size_t count) const {
            return getLHS().template packetPartial<AnyPacket>(index, count) - getRHS().template packetPartial<AnyPacket>(index, count);
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
    template<Vector V, class Executor>
    inline void VectorExpr<ExprType::Sub, T1, T2>::assign(V& v) const {
        constexpr bool FastAssign1 = Traits<T1>::FastAssign;
        constexpr bool FastAssign2 = Traits<T2>::FastAssign;
        if constexpr (FastAssign1) {
            getLHS().template assign<V, Executor>(v);
            v -= getRHS();
        }
        else if constexpr (FastAssign2) {
            (-getRHS()).template assign<V, Executor>(v);
            v += getLHS();
        }
        else
            Base::template assign<V, Executor>(v);
    }

    template<Vector T, Scalar U>
    [[nodiscard]] inline auto operator-(const T& v, const U& s) noexcept {
        return VectorExpr<ExprType::Sub, T, U>(v, s);
    }

    template<Scalar T, Vector U>
    [[nodiscard]] inline auto operator-(const T& v, const U& s) noexcept {
        return VectorExpr<ExprType::Sub, T, U>(v, s);
    }

    template<Vector T1, Vector T2>
    [[nodiscard]] inline auto operator-(const T1& v1, const T2& v2) noexcept {
        return VectorExpr<ExprType::Sub, T1, T2>(v1, v2);
    }
}
