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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.h"

namespace Physica {
    template<class U, class V>
    class VectorExpr<ExprID::Div, U, V>
            : public BinaryVectorExpr<ExprID::Div, U, V> {
        static_assert(Scalar<U> || Scalar<V>, "[Error]: Either type should be Scalar");

        using Base = BinaryVectorExpr<ExprID::Div, U, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        VectorExpr(U lhs, V rhs);
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept;
        [[nodiscard]] static CoDiff<T> operator()(const Scalar auto& lhs, std::random_access_iterator auto rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs, size_t count) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(const Scalar auto& lhs, std::random_access_iterator auto rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(const Scalar auto& lhs, std::random_access_iterator auto rhs, size_t count) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t s) const;
        [[nodiscard]] Tv calc_value(size_t s) const;

        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index, size_t count) const noexcept;
    };

    template<class U, class V>
    auto VectorExpr<ExprID::Div, U, V>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept -> CoDiff<T> {
        return *lhs / rhs;
    }

    template<class U, class V>
    auto VectorExpr<ExprID::Div, U, V>::operator()(const Scalar auto& lhs, std::random_access_iterator auto rhs) noexcept -> CoDiff<T> {
        assert(!(*rhs).isSubNormal() && "[Error]: Division overflow");
        return lhs / (*rhs);
    }

    template<class U, class V>
    template<int Size>
    auto VectorExpr<ExprID::Div, U, V>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs) noexcept -> SIMD<T, Size> {
        return lhs.template load<Size>() / SIMD<T, Size>(rhs);
    }

    template<class U, class V>
    template<int Size>
    auto VectorExpr<ExprID::Div, U, V>::operator()(std::random_access_iterator auto lhs, const Scalar auto& rhs, size_t count) noexcept -> SIMD<T, Size> {
        return lhs.template load<Size>(count) / SIMD<T, Size>(rhs);
    }

    template<class U, class V>
    template<int Size>
    auto VectorExpr<ExprID::Div, U, V>::operator()(const Scalar auto& lhs, std::random_access_iterator auto rhs) noexcept -> SIMD<T, Size> {
        auto div = rhs.template load<Size>();
        assert(!div.isSubNormal().horizontal_or() && "[Error]: Division overflow");
        return SIMD<T, Size>(lhs) / div;
    }

    template<class U, class V>
    template<int Size>
    auto VectorExpr<ExprID::Div, U, V>::operator()(const Scalar auto& lhs, std::random_access_iterator auto rhs, size_t count) noexcept -> SIMD<T, Size> {
        auto div = rhs.template load<Size>(count);
        for (size_t i = 0; i < count; ++i)
            assert(!div[i].isSubNormal() && "[Error]: Division overflow");
        return (SIMD<T, Size>(lhs) / div).cutoff(count);
    }

    template<class U, class V>
    VectorExpr<ExprID::Div, U, V>::VectorExpr(U lhs, V rhs) : Base(std::forward<U>(lhs), std::forward<V>(rhs)) {
        if constexpr (Vector<U>)
            assert(!Base::getRHS().isSubNormal() && "[Error]: Division overflow");
    }

    template<class U, class V>
    auto VectorExpr<ExprID::Div, U, V>::calc(size_t s) const -> CoDiff<T> {
        if constexpr (Vector<U>)
            return Base::getLHS().calc(s) / Base::getRHS();
        else
            return Base::getLHS() / Base::getRHS().calc(s);
    }

    template<class U, class V>
    auto VectorExpr<ExprID::Div, U, V>::calc_value(size_t s) const -> Tv {
        if constexpr (Vector<U>)
            return Base::getLHS().calc_value(s) / Base::getRHS().value();
        else
            return Base::getLHS().value() / Base::getRHS().calc_value(s);
    }

    template<class U, class V>
    template<int Size>
    auto VectorExpr<ExprID::Div, U, V>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        if constexpr (Vector<U>)
            return Base::getLHS().template packet<Size>(index) * SIMD<T, Size>(reciprocal(Base::getRHS()));
        else {
            auto div = Base::getRHS().template packet<Size>(index);
            assert(!div.isZero().horizontal_or() && "[Error]: Divide by zero");
            return SIMD<T, Size>(Base::getLHS()) / div;
        }
    }

    template<class U, class V>
    template<int Size>
    auto VectorExpr<ExprID::Div, U, V>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size> {
        if constexpr (Vector<U>)
            return Base::getLHS().template packet<Size>(index, count) * SIMD<T, Size>(reciprocal(Base::getRHS()));
        else {
            auto div = Base::getRHS().template packet<Size>(index, count);
            for (size_t i = 0; i < count; ++i)
                assert(!div[i].isZero() && "[Error]: Divide by zero");
            return (SIMD<T, Size>(Base::getLHS()) / div).cutoff(count);
        }
    }

    template<Vector V1, Vector V2>
    class VectorExpr<ExprID::Div, V1, V2>
            : public BinaryVectorExpr<ExprID::Div, V1, V2> {
        using Base = BinaryVectorExpr<ExprID::Div, V1, V2>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs, size_t count) noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t s) const;
        [[nodiscard]] Tv calc_value(size_t s) const;

        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index, size_t count) const noexcept;
    };

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Div, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept -> CoDiff<T> {
        assert(!(*rhs).isSubNormal() && "[Error]: Division overflow");
        return *lhs / *rhs;
    }

    template<Vector V1, Vector V2>
    template<int Size>
    auto VectorExpr<ExprID::Div, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs) noexcept -> SIMD<T, Size> {
        auto div = rhs.template load<Size>();
        assert(!div.isSubNormal().horizontal_or() && "[Error]: Division overflow");
        return lhs.template load<Size>() / div;
    }

    template<Vector V1, Vector V2>
    template<int Size>
    auto VectorExpr<ExprID::Div, V1, V2>::operator()(std::random_access_iterator auto lhs, std::random_access_iterator auto rhs, size_t count) noexcept -> SIMD<T, Size> {
        const auto pack1 = lhs.template load<Size>(count);
        const auto pack2 = rhs.template load<Size>(count);
        for (size_t i = 0; i < count; ++i)
            assert(!pack2[i].isSubNormal() && "[Error]: Division overflow");
        return (pack1 / pack2).cutoff(count);
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Div, V1, V2>::calc(size_t s) const -> CoDiff<T> {
        assert(!Base::getRHS().calc(s).isSubNormal() && "[Error]: Division overflow");
        return Base::getLHS().calc(s) / Base::getRHS().calc(s);
    }

    template<Vector V1, Vector V2>
    auto VectorExpr<ExprID::Div, V1, V2>::calc_value(size_t s) const -> Tv {
        assert(!Base::getRHS().calc_value(s).isSubNormal() && "[Error]: Division overflow");
        return Base::getLHS().calc_value(s) / Base::getRHS().calc_value(s);
    }

    template<Vector V1, Vector V2>
    template<int Size>
    auto VectorExpr<ExprID::Div, V1, V2>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        auto div = Base::getRHS().template packet<Size>(index);
        assert(!div.isZero().horizontal_or() && "[Error]: Divide by zero");
        return Base::getLHS().template packet<Size>(index) / div;
    }

    template<Vector V1, Vector V2>
    template<int Size>
    auto VectorExpr<ExprID::Div, V1, V2>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size> {
        const auto pack1 = Base::getLHS().template packet<Size>(index, count);
        const auto pack2 = Base::getRHS().template packet<Size>(index, count);
        for (size_t i = 0; i < count; ++i)
            assert(!pack2[i].isZero() && "[Error]: Divide by zero");
        return (pack1 / pack2).cutoff(count);
    }

    template<Vector V, Scalar U>
    [[nodiscard, gnu::always_inline]] auto operator/(V&& v, U&& x) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Div, V&&, U&&>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard, gnu::always_inline]] auto divide(U&& x, V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Div, U&&, V&&>(std::forward<U>(x), std::forward<V>(v));
    }

    template<Vector V1, Vector V2>
    [[nodiscard, gnu::always_inline]] auto divide(V1&& v1, V2&& v2) noexcept requires(!DeviceObj<V1> && !DeviceObj<V2>) {
        return VectorExpr<ExprID::Div, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
