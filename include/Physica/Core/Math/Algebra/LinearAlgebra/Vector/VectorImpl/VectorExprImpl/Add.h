/*
 * Copyright 2024-2026 Weibo He.
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
        [[nodiscard]] Pack packet(size_t index) const noexcept;
        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index, size_t count) const noexcept;

        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] auto values() const noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    private:
        consteval static bool lowerToFMA() noexcept;
    };

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Add, V, U>::calc(size_t s) const -> CoDiff<T> {
        if constexpr (lowerToFMA())
            return fma(getLHS().getLHS().calc(s), getLHS().getRHS().calc(s), getRHS());
        else
            return getLHS().calc(s) + getRHS();
    }

    template<Vector V, Scalar U>
    auto VectorExpr<ExprID::Add, V, U>::calc_value(size_t s) const -> Tv {
        if constexpr (lowerToFMA())
            return fma(getLHS().getLHS().calc_value(s), getLHS().getRHS().calc_value(s), getRHS().value());
        else
            return getLHS().calc_value(s) + getRHS().value();
    }

    template<Vector V, Scalar U>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Add, V, U>::packet(size_t index) const noexcept {
        if constexpr (lowerToFMA())
            return fma(getLHS().getLHS().template packet<Pack>(index),
                       getLHS().getRHS().template packet<Pack>(index),
                       Pack(getRHS()));
        else
            return getLHS().template packet<Pack>(index) + Pack(getRHS());
    }

    template<Vector V, Scalar U>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Add, V, U>::packet(size_t index, size_t count) const noexcept {
        if constexpr (lowerToFMA())
            return fma(getLHS().getLHS().template packet<Pack>(index, count),
                       getLHS().getRHS().template packet<Pack>(index, count),
                       Pack(getRHS(), count));
        else
            return getLHS().template packet<Pack>(index, count) + Pack(getRHS(), count);
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

    template<Vector V, Scalar U>
    consteval bool VectorExpr<ExprID::Add, V, U>::lowerToFMA() noexcept {
        if constexpr (instanceof_xt<VectorExpr, V>) {
            using V1 = std::remove_cvref_t<V>;
            using T1 = V1::ScalarType;
            using T2 = std::remove_cvref_t<U>;
            return (V1::getExprID() == ExprID::Mul) && std::same_as<T1, T2>;
        }
        return false;
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
        [[nodiscard]] Pack packet(size_t index, size_t count) const;

        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] auto values() const noexcept { return getLHS().values() + getRHS().values(); }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    private:
        consteval static bool lowerToFMA() noexcept;
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
            constexpr size_t SizeAtCompile = Base::maxSizeAtCompile(v);
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
        if constexpr (lowerToFMA()) {
            if constexpr (Scalar<decltype(getLHS().getRHS())>)
                return fma(getLHS().getLHS().calc(s), getLHS().getRHS(), getRHS().calc(s));
            else
                return fma(getLHS().getLHS().calc(s), getLHS().getRHS().calc(s), getRHS().calc(s));
        }
        else
            return getLHS().calc(s) + getRHS().calc(s);
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Add, V1, V2>::calc_value(size_t s) const -> Tv {
        if constexpr (lowerToFMA()) {
            if constexpr (Scalar<decltype(getLHS().getRHS())>)
                return fma(getLHS().getLHS().calc(s), getLHS().getRHS().value(), getRHS().calc(s));
            else
                return fma(getLHS().getLHS().calc_value(s), getLHS().getRHS().calc_value(s), getRHS().calc_value(s));
        }
        else
            return getLHS().calc_value(s) + getRHS().calc_value(s);
    }

    template<Vector V1, Vector V2>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Add, V1, V2>::packet(size_t index) const {
        if constexpr (lowerToFMA()) {
            if constexpr (Scalar<decltype(getLHS().getRHS())>) {
                return fma(getLHS().getLHS().template packet<Pack>(index),
                           Pack(getLHS().getRHS()),
                           getRHS().template packet<Pack>(index));
            }
            else {
                return fma(getLHS().getLHS().template packet<Pack>(index),
                           getLHS().getRHS().template packet<Pack>(index),
                           getRHS().template packet<Pack>(index));
            }
        }
        else
            return getLHS().template packet<Pack>(index) + getRHS().template packet<Pack>(index);
    }

    template<Vector V1, Vector V2>
    template<Packet Pack>
    Pack VectorExpr<ExprID::Add, V1, V2>::packet(size_t index, size_t count) const {
        if constexpr (lowerToFMA()) {
            if constexpr (Scalar<decltype(getLHS().getRHS())>) {
                return fma(getLHS().getLHS().template packet<Pack>(index, count),
                           Pack(getLHS().getRHS()).cutoff(count),
                           getRHS().template packet<Pack>(index, count));
            }
            else {
                return fma(getLHS().getLHS().template packet<Pack>(index, count),
                           getLHS().getRHS().template packet<Pack>(index, count),
                           getRHS().template packet<Pack>(index, count));
            }
        }
        else
            return getLHS().template packet<Pack>(index, count) + getRHS().template packet<Pack>(index, count);
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

    template<Vector V1, Vector V2>
    consteval bool VectorExpr<ExprID::Add, V1, V2>::lowerToFMA() noexcept {
        if constexpr (instanceof_xt<VectorExpr, V1>) {
            using T1 = std::remove_cvref_t<V1>::ScalarType;
            using T2 = std::remove_cvref_t<V2>::ScalarType;
            return (std::remove_cvref_t<V1>::getExprID() == ExprID::Mul) && std::same_as<T1, T2>;
        }
        return false;
    }

    template<Vector V, Scalar U>
    [[nodiscard]] auto operator+(V&& v, U&& x) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Add, V&&, U&&>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard]] auto operator+(U&& x, V&& v) noexcept requires(!DeviceObj<V>) {
        return std::forward<V>(v) + std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto operator+(V1&& v1, V2&& v2) noexcept requires(!DeviceObj<V1> && !DeviceObj<V2>) {
        return VectorExpr<ExprID::Add, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}

#ifdef PHYSICA_MKL
    #include "MKL/Add.h"
#endif
