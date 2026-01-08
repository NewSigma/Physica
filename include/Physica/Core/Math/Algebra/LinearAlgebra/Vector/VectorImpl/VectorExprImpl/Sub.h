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
    template<class T1, class T2>
    class VectorExpr<ExprID::Sub, T1, T2>
            : public BinaryVectorExpr<ExprID::Sub, T1, T2> {
        static_assert(Scalar<T1> || Scalar<T2>, "[Error]: Either type should be Scalar");

        using Base = BinaryVectorExpr<ExprID::Sub, T1, T2>;
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

        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] auto values() const noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<class T1, class T2>
    auto VectorExpr<ExprID::Sub, T1, T2>::calc(size_t s) const -> CoDiff<T> {
        if constexpr (Vector<T1>)
            return getLHS().calc(s) - getRHS();
        else
            return getLHS() - getRHS().calc(s);
    }

    template<class T1, class T2>
    auto VectorExpr<ExprID::Sub, T1, T2>::calc_value(size_t s) const -> Tv {
        if constexpr (Vector<T1>)
            return getLHS().calc_value(s) - getRHS().value();
        else
            return getLHS().value() - getRHS().calc_value(s);
    }

    template<class T1, class T2>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Sub, T1, T2>::packet(size_t index) const {
        if constexpr (Vector<T1>)
            return getLHS().template packet<Pack>(index) - Pack(getRHS());
        else
            return Pack(getLHS()) - getRHS().template packet<Pack>(index);
    }

    template<class T1, class T2>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Sub, T1, T2>::packetPartial(size_t index, size_t count) const {
        if constexpr (Vector<T1>)
            return getLHS().template packetPartial<Pack>(index, count) - Pack(getRHS(), count);
        else
            return Pack(getLHS(), count) - getRHS().template packetPartial<Pack>(index, count);
    }

    template<class T1, class T2>
    void VectorExpr<ExprID::Sub, T1, T2>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        if constexpr (Scalar<decltype(grad)>) {
            const auto& g = grad.value();
            if constexpr (Vector<T1>) {
                if constexpr (ReverseDiff<T1>)
                    lhs.reverse(g);
                if constexpr (ReverseDiff<T2>)
                    rhs.reverse(-g * Tv(Base::getLength()));
            }
            else {
                if constexpr (ReverseDiff<T1>)
                    lhs.reverse(g * Tv(Base::getLength()));
                if constexpr (ReverseDiff<T2>)
                    rhs.reverse(-g);
            }
        }
        else {
            static_assert(Vector<decltype(grad)>, "[Error]: Unexpected type");
            const auto& g = grad.values();
            if constexpr (Vector<T1>) {
                if constexpr (ReverseDiff<T1>)
                    lhs.reverse(g);
                if constexpr (ReverseDiff<T2>)
                    rhs.reverse(-g.sum());
            }
            else {
                if constexpr (ReverseDiff<T1>)
                    lhs.reverse(g.sum());
                if constexpr (ReverseDiff<T2>)
                    rhs.reverse(-g);
            }
        }
    }

    template<class T1, class T2>
    auto VectorExpr<ExprID::Sub, T1, T2>::values() const noexcept {
        if constexpr (Vector<T1>)
            return getLHS().values() - getRHS().value();
        else
            return getLHS().value() - getRHS().values();
    }

    template<Vector V1, Vector V2>
    class VectorExpr<ExprID::Sub, V1, V2>
            : public BinaryVectorExpr<ExprID::Sub, V1, V2> {
        using Base = BinaryVectorExpr<ExprID::Sub, V1, V2>;
    public:
        using Base::isComplex;
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

        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] auto values() const noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector V1, Vector V2>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Sub, V1, V2>::assign(Vector auto&& v) const {
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
        else {
            using V = std::remove_cvref<decltype(v)>::type;
            constexpr size_t SizeAtCompile = Base::maxSizeAtCompile(v);
            constexpr size_t Critical = HostDevAttr::LineSizeL1D / sizeof(T);
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
    auto VectorExpr<ExprID::Sub, V1, V2>::calc(size_t s) const -> CoDiff<T> {
        return getLHS().calc(s) - getRHS().calc(s);
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Sub, V1, V2>::calc_value(size_t s) const -> Tv {
        return getLHS().calc_value(s) - getRHS().calc_value(s);
    }

    template<Vector V1, Vector V2>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Sub, V1, V2>::packet(size_t index) const {
        return getLHS().template packet<Pack>(index) - getRHS().template packet<Pack>(index);
    }

    template<Vector V1, Vector V2>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Sub, V1, V2>::packetPartial(size_t index, size_t count) const {
        return getLHS().template packetPartial<Pack>(index, count) - getRHS().template packetPartial<Pack>(index, count);
    }

    template<Vector V1, Vector V2>
    void VectorExpr<ExprID::Sub, V1, V2>::reverse(const Vector auto& grad) const noexcept {
        static_assert(Base::isReverseDiff);
        const auto& g = grad.values();
        if constexpr (ReverseDiff<V1>)
            Base::getLHS().reverse(g);
        if constexpr (ReverseDiff<V2>)
            Base::getRHS().reverse(g);
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Sub, V1, V2>::values() const noexcept {
        return getLHS().values() - getRHS().values();
    }

    template<Vector V, Scalar U>
    [[nodiscard]] auto operator-(V&& v, U&& x) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprID::Sub, V&&, U&&>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard]] auto operator-(U&& x, V&& v) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprID::Sub, U&&, V&&>(std::forward<U>(x), std::forward<V>(v));
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto operator-(V1&& v1, V2&& v2) noexcept requires(!CUDA<V1> && !CUDA<V2>) {
        return VectorExpr<ExprID::Sub, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}

#ifdef PHYSICA_MKL
    #include "MKL/Sub.h"
#endif
