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
    template<Vector V>
    class VectorExpr<ExprID::Abs, V> : public UnitaryVectorExpr<ExprID::Abs, V> {
        using This = VectorExpr<ExprID::Abs, V>;
        using Base = UnitaryVectorExpr<ExprID::Abs, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        constexpr static bool isComplexV = std::remove_cvref_t<V>::isComplex;
        constexpr static bool isReverseDiff = Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] static CoDiff<T> operator()(std::random_access_iterator auto input) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input) noexcept;
        template<int Size>
        [[nodiscard]] static SIMD<T, Size> operator()(std::random_access_iterator auto input, size_t count) noexcept;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const;
        void assign_mkl(Vector auto& v) const noexcept;

        [[nodiscard]] CoDiff<T> calc(size_t index) const { return abs(Base::getExpr().calc(index)); }
        [[nodiscard]] Tv calc_value(size_t index) const { return abs(Base::getExpr().calc_value(index)); }
        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] T max() const;

        [[nodiscard]] auto values(this auto&&) noexcept;
    };

    template<Vector V>
    auto VectorExpr<ExprID::Abs, V>::operator()(std::random_access_iterator auto input) noexcept -> CoDiff<T> {
        return abs(*input);
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Abs, V>::operator()(std::random_access_iterator auto input) noexcept -> SIMD<T, Size> {
        if constexpr (isComplexV)
            return sqrt(SquaredNormVector<V>::template operator()<Size>(input));
        else
            return abs(input.template load<Size>());
    }

    template<Vector V>
    template<int Size>
    auto VectorExpr<ExprID::Abs, V>::operator()(std::random_access_iterator auto input, size_t count) noexcept -> SIMD<T, Size> {
        if constexpr (isComplexV)
            return sqrt(SquaredNormVector<V>::template operator()<Size>(input, count));
        else
            return abs(input.template load<Size>(count));
    }

    template<Vector V>
    template<ExecutePolicy P>
    void VectorExpr<ExprID::Abs, V>::assign(Vector auto&& v) const {
        Base::template assign_base<P>(v);
    }

    template<Vector V>
    void VectorExpr<ExprID::Abs, V>::reverse(const auto& grad) const noexcept {
        static_assert(isReverseDiff);
        using U = decltype(grad);
        const auto& expr = Base::getExpr();
        if constexpr (Scalar<U>)
            expr.reverse(unit(expr.values()) * grad);
        else {
            static_assert(Vector<U>, "[Error]: Unexpected type");
            expr.reverse(hadamard(unit(expr.values()), grad));
        }
    }

    template<Vector V>
    auto VectorExpr<ExprID::Abs, V>::max() const -> T {
        if constexpr (isComplexV)
            return sqrt(Base::getExpr().squaredNorms().max());
        else
            return Base::max();
    }

    template<Vector V>
    auto VectorExpr<ExprID::Abs, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return abs(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] auto abs(V&& v) noexcept requires(!DeviceObj<V>) {
        return VectorExpr<ExprID::Abs, V&&>(std::forward<V>(v));
    }
}

#ifdef PHYSICA_MKL
    #include "MKL/Abs.h"
#endif
