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
    template<Vector V, Scalar U>
    class VectorExpr<ExprType::Add, V, U>
            : public BinaryVectorExpr<ExprType::Add, V, U> {
        using Base = BinaryVectorExpr<ExprType::Add, V, U>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t s) const {
            return Base::getLHS().calc(s) + Base::getRHS();
        }

        [[nodiscard]] Tv calc_value(size_t s) const {
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

        void reverse(const Vector auto& grad_) const noexcept requires(isReverseDiff) {
            const auto& grad = grad_.values();
            if constexpr (ReverseDiff<V>)
                Base::getLHS().reverse(grad);
            if constexpr (ReverseDiff<U>)
                Base::getRHS().reverse(grad.sum());
        }
    };

    template<Vector V1, Vector V2>
    class VectorExpr<ExprType::Add, V1, V2>
            : public BinaryVectorExpr<ExprType::Add, V1, V2> {
        using Base = BinaryVectorExpr<ExprType::Add, V1, V2>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& v) const;

        [[nodiscard]] CoDiff<T> calc(size_t s) const {
            return getLHS().calc(s) + getRHS().calc(s);
        }

        [[nodiscard]] Tv calc_value(size_t s) const {
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

        void reverse(const auto& grad_) const noexcept requires(isReverseDiff);

        auto values() const noexcept { return Base::getLHS().values() + Base::getRHS().values(); }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector V1, Vector V2>
    template<ExecutePolicy P>
    void VectorExpr<ExprType::Add, V1, V2>::assign(Vector auto& v) const {
        constexpr bool FastAssign1 = Traits<std::remove_cvref_t<V1>>::FastAssign;
        constexpr bool FastAssign2 = Traits<std::remove_cvref_t<V2>>::FastAssign;
        if constexpr (FastAssign1) {
            getLHS().template assign<P>(v);
            v += getRHS();
        }
        else if constexpr (FastAssign2) {
                getRHS().template assign<P>(v);
                v += getLHS();
        }
        else
            Base::template assign<P>(v);
    }

    template<Vector V1, Vector V2>
    void VectorExpr<ExprType::Add, V1, V2>::reverse(const auto& grad_) const noexcept requires(isReverseDiff) {
        using U = decltype(grad_);
        if constexpr (Scalar<U>) {
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(grad_);
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(grad_);
        }
        else {
            static_assert(Vector<U>);
            const auto& grad = grad_.values();
            assert(grad.getLength() == Base::getLength());
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(grad);
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(grad);
        }
    }

    template<Vector V, Scalar U>
    [[nodiscard]] auto operator+(V&& v, U&& x) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Add, V&&, U&&>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard]] auto operator+(U&& x, V&& v) noexcept requires(!CUDA<V>) {
        return std::forward<V>(v) + std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto operator+(V1&& v1, V2&& v2) noexcept requires(!CUDA<V1> && !CUDA<V2>) {
        return VectorExpr<ExprType::Add, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
