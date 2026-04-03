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
    template<Matrix M, Scalar U>
    class MatrixExpr<ExprID::Add, M, U>
            : public BinaryMatrixExpr<ExprID::Add, M, U> {
        using Base = BinaryMatrixExpr<ExprID::Add, M, U>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] auto values(this auto&& self) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M, Scalar U>
    void MatrixExpr<ExprID::Add, M, U>::assign(Matrix auto& target) const {
        if constexpr (MatrixMajor::isSameMajor<M, decltype(target)>())
            (getLHS().flatten() + getRHS()).assign(target.flatten());
        else
            Base::assign(target);
    }

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprID::Add, M, U>::calc(size_t row, size_t col) const -> T {
        return getLHS().calc(row, col) + getRHS();
    }

    template<Matrix M, Scalar U>
    void MatrixExpr<ExprID::Add, M, U>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(Base::isReverseDiff());
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        if constexpr (Diffable<M>)
            lhs.reverse(grad);
        if constexpr (Diffable<U>)
            rhs.reverse(grad.sum());
    }

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprID::Add, M, U>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() + std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Vector V>
    class MatrixExpr<ExprID::Add, M, V>
            : public BinaryMatrixExpr<ExprID::Add, M, V> {
        using Base = BinaryMatrixExpr<ExprID::Add, M, V>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;

        [[nodiscard]] auto values(this auto&& self) noexcept;
    };

    template<Matrix M, Vector V>
    auto MatrixExpr<ExprID::Add, M, V>::calc(size_t row, size_t col) const -> T {
        return Base::getLHS().calc(row, col) + Base::getRHS().calc(row);
    }

    template<Matrix M, Vector V>
    auto MatrixExpr<ExprID::Add, M, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() + std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M1, Matrix M2>
    class MatrixExpr<ExprID::Add, M1, M2>
            : public BinaryMatrixExpr<ExprID::Add, M1, M2> {
        using Base = BinaryMatrixExpr<ExprID::Add, M1, M2>;
        using This = MatrixExpr<ExprID::Add, M1, M2>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] auto values(this auto&& self) noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprID::Add, M1, M2>::assign(Matrix auto& target) const {
        constexpr bool SameMajor1 = MatrixMajor::isSameMajor<M1, decltype(target)>();
        constexpr bool SameMajor2 = MatrixMajor::isSameMajor<M2, decltype(target)>();
        if constexpr (SameMajor1 && SameMajor2)
            (getLHS().flatten() + getRHS().flatten()).assign(target.flatten());
        else
            Base::assign(target);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprID::Add, M1, M2>::calc(size_t row, size_t col) const -> T {
        return Base::getLHS().calc(row, col) + Base::getRHS().calc(row, col);
    }

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprID::Add, M1, M2>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(Base::isReverseDiff());
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        if constexpr (Diffable<M1>)
            lhs.reverse(grad);
        if constexpr (Diffable<M2>)
            rhs.reverse(grad);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprID::Add, M1, M2>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() + std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] auto operator+(M&& m, U&& x) noexcept requires(!DeviceObj<M>) {
        return MatrixExpr<ExprID::Add, M&&, U&&>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] auto operator+(U&& x, M&& m) noexcept requires(!DeviceObj<M>) {
        return std::forward<M>(m) + std::forward<U>(x);
    }

    template<Matrix M, Vector V>
    [[nodiscard, gnu::always_inline]] auto operator+(M&& m, V&& x) noexcept requires(!DeviceObj<M> && !DeviceObj<V>) {
        return MatrixExpr<ExprID::Add, M&&, V&&>(std::forward<M>(m), std::forward<V>(x));
    }

    template<Matrix M, Vector V>
    [[nodiscard, gnu::always_inline]] auto operator+(V&& x, M&& m) noexcept requires(!DeviceObj<M> && !DeviceObj<V>) {
        return std::forward<M>(m) + std::forward<V>(x);
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard, gnu::always_inline]] auto operator+(M1&& m1, M2&& m2) noexcept requires(!DeviceObj<M1> && !DeviceObj<M2>) {
        if constexpr (!canonicalized(m1, m2))
            return std::forward<M2>(m2) + std::forward<M1>(m1);
        else
            return MatrixExpr<ExprID::Add, M1&&, M2&&>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}
