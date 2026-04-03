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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.h"

namespace Physica {
    template<Matrix M1, Matrix M2>
    class MatrixExpr<ExprID::Sub, M1, M2>
            : public BinaryMatrixExpr<ExprID::Sub, M1, M2> {
        using Base = BinaryMatrixExpr<ExprID::Sub, M1, M2>;
        using This = MatrixExpr<ExprID::Sub, M1, M2>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        void assign(Matrix auto&& target) const;

        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;

        using Base::reverse;
        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprID::Sub, M1, M2>::assign(Matrix auto&& target) const {
        constexpr bool SameMajor1 = MatrixMajor::isSameMajor<M1, decltype(target)>();
        constexpr bool SameMajor2 = MatrixMajor::isSameMajor<M2, decltype(target)>();
        if constexpr (SameMajor1 && SameMajor2)
            (getLHS().flatten() - getRHS().flatten()).assign(target.flatten());
        else
            Base::assign(target);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprID::Sub, M1, M2>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return getLHS().calc(row, col) - getRHS().calc(row, col);
    }

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprID::Sub, M1, M2>::reverse(const auto& grad) const noexcept {
        if constexpr (ReverseDiff<M1>)
            getLHS().reverse(grad);
        if constexpr (ReverseDiff<M2>)
            getRHS().reverse(-grad);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprID::Sub, M1, M2>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() - std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] auto operator-(M&& m, U&& x) noexcept requires(!DeviceObj<M>) {
        return std::forward<M>(m) + (-std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] auto operator-(U&& x, M&& m) noexcept requires(!DeviceObj<M>) {
        return std::forward<U>(x) + (-std::forward<M>(m));
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard, gnu::always_inline]] auto operator-(M1&& m1, M2&& m2) noexcept requires(!DeviceObj<M1> && !DeviceObj<M2>) {
        return MatrixExpr<ExprID::Sub, M1&&, M2&&>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}
