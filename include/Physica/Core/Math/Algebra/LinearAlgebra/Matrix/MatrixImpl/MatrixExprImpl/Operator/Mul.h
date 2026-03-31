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
    class MatrixExpr<ExprID::Mul, M, U>
            : public BinaryMatrixExpr<ExprID::Mul, M, U> {
        using Base = BinaryMatrixExpr<ExprID::Mul, M, U>;
        using This = MatrixExpr<ExprID::Mul, M, U>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        using Base::operator*;
        [[nodiscard]] auto operator*(Scalar auto x) const noexcept;
        [[nodiscard]] auto operator-(this auto&&) noexcept;
        /* Operations */
        void assign(Matrix auto&& target) const;

        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        [[nodiscard]] T trace() const { return getLHS().trace() * getRHS(); }
        [[nodiscard]] T sum() const { return getLHS().sum() * getRHS(); }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprID::Mul, M, U>::operator*(Scalar auto x) const noexcept {
        return getLHS() * (getRHS() * x);
    }

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprID::Mul, M, U>::operator-(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS() * (-std::forward<Self>(self).getRHS());
    }

    template<Matrix M, Scalar U>
    void MatrixExpr<ExprID::Mul, M, U>::assign(Matrix auto&& target) const {
        if constexpr (MatrixMajor::isSameMajor<M, decltype(target)>())
            (getLHS().flatten() * getRHS()).assign(target.flatten());
        else
            Base::assign(target);
    }

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprID::Mul, M, U>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return getLHS().calc(row, col) * getRHS();
    }

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprID::Mul, M, U>::calc_value(size_t row, size_t col) const -> Tv {
        return getLHS().calc_value(row, col) * getRHS().value();
    }

    template<Matrix M1, Matrix M2>
    class MatrixExpr<ExprID::Mul, M1, M2>
            : public BinaryMatrixExpr<ExprID::Mul, M1, M2> {
        using Base = BinaryMatrixExpr<ExprID::Mul, M1, M2>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        void assign(Matrix auto&& target) const;

        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprID::Mul, M1, M2>::assign(Matrix auto&& target) const {
        constexpr bool SameMajor1 = MatrixMajor::isSameMajor<M1, decltype(target)>();
        constexpr bool SameMajor2 = MatrixMajor::isSameMajor<M2, decltype(target)>();
        if constexpr (SameMajor1 && SameMajor2)
            hadamard(getLHS().flatten(), getRHS().flatten()).assign(target.flatten());
        else
            Base::assign(target);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprID::Mul, M1, M2>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return getLHS().calc(row, col) * getRHS().calc(row, col);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprID::Mul, M1, M2>::calc_value(size_t row, size_t col) const -> Tv {
        return getLHS().calc_value(row, col) * getRHS().calc_value(row, col);
    }

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprID::Mul, M1, M2>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        if constexpr (Diffable<M1>)
            lhs.reverse(hadamard(rhs.values(), grad));
        if constexpr (Diffable<M2>)
            rhs.reverse(hadamard(lhs.values(), grad));
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] auto operator*(M&& m, U&& x) noexcept requires(!DeviceObj<M>) {
        return MatrixExpr<ExprID::Mul, M&&, U&&>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] auto operator*(U&& x, M&& m) noexcept requires(!DeviceObj<M>) {
        return std::forward<M>(m) * std::forward<U>(x);
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard, gnu::always_inline]] auto hadamard(M1&& m1, M2&& m2) noexcept requires(!DeviceObj<M1> && !DeviceObj<M2>) {
        if constexpr (!canonicalized(m1, m2))
            return hadamard(std::forward<M2>(m2), std::forward<M1>(m1));
        else
            return MatrixExpr<ExprID::Mul, M1&&, M2&&>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}

#include "GEVM/Mul.h"
