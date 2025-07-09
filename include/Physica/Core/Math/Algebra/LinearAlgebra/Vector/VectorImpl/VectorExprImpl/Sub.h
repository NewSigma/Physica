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
    template<class T, class U>
    class VectorExpr<ExprType::Sub, T, U>
            : public BinaryVectorExpr<ExprType::Sub, T, U> {
        static_assert(Scalar<T> || Scalar<U>, "[Error]: Either type should be Scalar");

        using Base = BinaryVectorExpr<ExprType::Sub, T, U>;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
    protected:
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const {
            if constexpr (Vector<T>)
                return Base::getLHS().calc(s) - Base::getRHS();
            else
                return Base::getLHS() - Base::getRHS().calc(s);
        }

        [[nodiscard]] Tv calc_value(size_t s) const {
            if constexpr (Vector<T>)
                return Base::getLHS().calc_value(s) - Base::getRHS().value();
            else
                return Base::getLHS().value() - Base::getRHS().calc_value(s);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            if constexpr (Vector<T>)
                return Base::getLHS().template packet<Pack>(index) - Pack(Base::getRHS());
            else
                return Pack(Base::getLHS()) - Base::getRHS().template packet<Pack>(index);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            if constexpr (Vector<T>)
                return Base::getLHS().template packetPartial<Pack>(index, count) - Pack(Base::getRHS(), count);
            else
                return Pack(Base::getLHS(), count) - Base::getRHS().template packetPartial<Pack>(index, count);
        }

        void reverse(const auto& grad) const noexcept requires(isReverseDiff);
    };

    template<class T, class U>
    void VectorExpr<ExprType::Sub, T, U>::reverse(const auto& grad) const noexcept requires(isReverseDiff) {
        const auto& lhs = Base::getLHS();
        const auto& rhs = Base::getRHS();
        if constexpr (Scalar<decltype(grad)>) {
            const auto& g = grad.value();
            if constexpr (Vector<T>) {
                if constexpr (ReverseDiff<T>)
                    lhs.reverse(g);
                if constexpr (ReverseDiff<U>)
                    rhs.reverse(-g * Tv(Base::getLength()));
            }
            else {
                if constexpr (ReverseDiff<T>)
                    lhs.reverse(g * Tv(Base::getLength()));
                if constexpr (ReverseDiff<U>)
                    rhs.reverse(-g);
            }
        }
        else {
            static_assert(Vector<decltype(grad)>, "[Error]: Unexpected type");
            const auto& g = grad.values();
            if constexpr (Vector<T>) {
                if constexpr (ReverseDiff<T>)
                    lhs.reverse(g);
                if constexpr (ReverseDiff<U>)
                    rhs.reverse(-g.sum());
            }
            else {
                if constexpr (ReverseDiff<T>)
                    lhs.reverse(g.sum());
                if constexpr (ReverseDiff<U>)
                    rhs.reverse(-g);
            }
        }
    }

    template<Vector V1, Vector V2>
    class VectorExpr<ExprType::Sub, V1, V2>
            : public BinaryVectorExpr<ExprType::Sub, V1, V2> {
        using Base = BinaryVectorExpr<ExprType::Sub, V1, V2>;
    public:
        using typename Base::ScalarType;
        using Base::isReverseDiff;
    protected:
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& v) const;

        [[nodiscard]] CoDiff<ScalarType> calc(size_t s) const {
            return getLHS().calc(s) - getRHS().calc(s);
        }

        [[nodiscard]] Tv calc_value(size_t s) const {
            return getLHS().calc_value(s) - getRHS().calc_value(s);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return getLHS().template packet<Pack>(index) - getRHS().template packet<Pack>(index);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return getLHS().template packetPartial<Pack>(index, count) - getRHS().template packetPartial<Pack>(index, count);
        }

        void reverse(const Vector auto& grad) const noexcept requires(isReverseDiff) {
            const auto& g = grad.values();
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(g);
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(g);
        }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector V1, Vector V2>
    template<ExecutePolicy P>
    void VectorExpr<ExprType::Sub, V1, V2>::assign(Vector auto& v) const {
        constexpr bool FastAssign1 = Traits<std::remove_cvref_t<V1>>::FastAssign;
        constexpr bool FastAssign2 = Traits<std::remove_cvref_t<V2>>::FastAssign;
        if constexpr (FastAssign1) {
            getLHS().template assign<P>(v);
            v -= getRHS();
        }
        else if constexpr (FastAssign2) {
            static_assert(Traits<decltype(-getRHS())>::FastAssign, "[Debug]: Fast minus implementation is missing");
            (-getRHS()).template assign<P>(v);
            v += getLHS();
        }
        else
            Base::template assign<P>(v);
    }

    template<Vector T, Scalar U>
    [[nodiscard]] auto operator-(T&& v, U&& x) noexcept requires(!CUDA<T>) {
        return VectorExpr<ExprType::Sub, T&&, U&&>(std::forward<T>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector T>
    [[nodiscard]] auto operator-(U&& x, T&& v) noexcept requires(!CUDA<T>) {
        return VectorExpr<ExprType::Sub, U&&, T&&>(std::forward<U>(x), std::forward<T>(v));
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto operator-(V1&& v1, V2&& v2) noexcept requires(!CUDA<V1> && !CUDA<V2>) {
        return VectorExpr<ExprType::Sub, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
