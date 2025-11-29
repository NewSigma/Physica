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
    class VectorExpr<ExprID::Add, V, U>
            : public BinaryVectorExpr<ExprID::Add, V, U> {
        using Base = BinaryVectorExpr<ExprID::Add, V, U>;
    public:
        using Base::isReverseDiff;
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

        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] auto values() const noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Add, V, U>::calc(size_t s) const -> CoDiff<T> {
        return Base::getLHS().calc(s) + Base::getRHS();
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Add, V, U>::calc_value(size_t s) const -> Tv {
        return Base::getLHS().calc_value(s) + Base::getRHS().value();
    }

    template<Vector V, Scalar U>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Add, V, U>::packet(size_t index) const {
        return Base::getLHS().template packet<Pack>(index) + Pack(Base::getRHS());
    }

    template<Vector V, Scalar U>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Add, V, U>::packetPartial(size_t index, size_t count) const {
        return Base::getLHS().template packetPartial<Pack>(index, count) + Pack(Base::getRHS(), count);
    }

    template<Vector V, Scalar U>
    void VectorExpr<ExprID::Add, V, U>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& g = grad.values();
        if constexpr (ReverseDiff<V>)
            Base::getLHS().reverse(g);
        if constexpr (ReverseDiff<U>)
            Base::getRHS().reverse(g.sum());
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Add, V, U>::values() const noexcept {
        return getLHS().values() + getRHS().value();
    }

    template<Vector V1, Vector V2>
    class VectorExpr<ExprID::Add, V1, V2>
            : public BinaryVectorExpr<ExprID::Add, V1, V2> {
        using Base = BinaryVectorExpr<ExprID::Add, V1, V2>;
    public:
        using Base::isComplex;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tc;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const;
        void assign_mkl(Vector auto& v) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t s) const;
        [[nodiscard]] Tv calc_value(size_t s) const;

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;

        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] auto values() const noexcept { return getLHS().values() + getRHS().values(); }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector V1, Vector V2>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Add, V1, V2>::assign(Vector auto&& v) const {
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
        else {
            using V = std::remove_cvref<decltype(v)>::type;
            constexpr size_t SizeAtCompile = Base::template maxSizeAtCompile<decltype(v)>();
            constexpr size_t Critical = 256;
            constexpr bool UseMKL1 = Internal::EnableMKL<V1, V>::value;
            constexpr bool UseMKL2 = Internal::EnableMKL<V2, V>::value;
            constexpr bool UseMKL3 = SizeAtCompile == Dynamic || SizeAtCompile > Critical;
            constexpr bool UseMKL = UseMKL1 && UseMKL2 && UseMKL3 && (T::Prec != Float64);
            if constexpr (UseMKL) {
                if (Base::getLength() > Critical)
                    assign_mkl(v);
                else
                    Base::template assign_base<P>(v);
            }
            else
                Base::template assign_base<P>(v);
        }
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Add, V1, V2>::calc(size_t s) const -> CoDiff<T> {
        return getLHS().calc(s) + getRHS().calc(s);
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Add, V1, V2>::calc_value(size_t s) const -> Tv {
        return Base::getLHS().calc_value(s) + Base::getRHS().calc_value(s);
    }

    template<Vector V1, Vector V2>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Add, V1, V2>::packet(size_t index) const {
        return getLHS().template packet<Pack>(index) + getRHS().template packet<Pack>(index);
    }

    template<Vector V1, Vector V2>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Add, V1, V2>::packetPartial(size_t index, size_t count) const {
        return getLHS().template packetPartial<Pack>(index, count) + getRHS().template packetPartial<Pack>(index, count);
    }

    template<Vector V1, Vector V2>
    void VectorExpr<ExprID::Add, V1, V2>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff);
        using U = decltype(grad);
        if constexpr (Scalar<U>) {
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(grad);
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(grad);
        }
        else {
            static_assert(Vector<U>);
            const auto& g = grad.values();
            assert(g.getLength() == Base::getLength());
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(g);
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(g);
        }
    }

    template<Vector V, Scalar U>
    [[nodiscard]] auto operator+(V&& v, U&& x) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprID::Add, V&&, U&&>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard]] auto operator+(U&& x, V&& v) noexcept requires(!CUDA<V>) {
        return std::forward<V>(v) + std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto operator+(V1&& v1, V2&& v2) noexcept requires(!CUDA<V1> && !CUDA<V2>) {
        return VectorExpr<ExprID::Add, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}

#ifdef PHYSICA_MKL
    #include "MKL/Add.h"
#endif
