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

#include "../MatrixExpr.h"

namespace Physica {
    template<Matrix M, Scalar U>
    class MatrixExpr<ExprType::Mul, M, U>
            : public BinaryMatrixExpr<ExprType::Mul, M, U> {
        using Base = BinaryMatrixExpr<ExprType::Mul, M, U>;
        using This = MatrixExpr<ExprType::Mul, M, U>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] auto operator-() const& noexcept;
        [[nodiscard]] auto operator-() && noexcept;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        [[nodiscard]] T trace() const { return getLHS().trace() * getRHS(); }
        [[nodiscard]] T sum() const { return getLHS().sum() * getRHS(); }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprType::Mul, M, U>::operator-() const& noexcept {
        return getLHS() * (-getRHS());
    }

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprType::Mul, M, U>::operator-() && noexcept {
        return std::move(getLHS()) * (-getRHS());
    }

    template<Matrix M, Scalar U>
    void MatrixExpr<ExprType::Mul, M, U>::assign(Matrix auto& target) const {
        if constexpr (MatrixOption::isSameMajor<M, decltype(target)>())
            (getLHS().flatten() * getRHS()).assign(target.flatten());
        else
            Base::assign(target);
    }

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprType::Mul, M, U>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return getLHS().calc(row, col) * getRHS();
    }

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprType::Mul, M, U>::calc_value(size_t row, size_t col) const -> Tv {
        return getLHS().calc_value(row, col) * getRHS().value();
    }

    template<Matrix M1, Matrix M2>
    class MatrixExpr<ExprType::Mul, M1, M2>
            : public BinaryMatrixExpr<ExprType::Mul, M1, M2> {
        using Base = BinaryMatrixExpr<ExprType::Mul, M1, M2>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        void assign(Matrix auto& target) const;

        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprType::Mul, M1, M2>::assign(Matrix auto& target) const {
        constexpr bool SameMajor1 = MatrixOption::isSameMajor<M1, decltype(target)>();
        constexpr bool SameMajor2 = MatrixOption::isSameMajor<M2, decltype(target)>();
        if constexpr (SameMajor1 && SameMajor2)
            hadamard(getLHS().flatten(), getRHS().flatten()).assign(target.flatten());
        else
            Base::assign(target);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprType::Mul, M1, M2>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return getLHS().calc(row, col) * getRHS().calc(row, col);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprType::Mul, M1, M2>::calc_value(size_t row, size_t col) const -> Tv {
        return getLHS().calc_value(row, col) * getRHS().calc_value(row, col);
    }

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprType::Mul, M1, M2>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        if constexpr (Diffable<M1>)
            lhs.reverse(hadamard(rhs.values(), grad));
        if constexpr (Diffable<M2>)
            rhs.reverse(hadamard(lhs.values(), grad));
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] auto operator*(M&& m, U&& x) noexcept requires(!CUDA<M>) {
        return MatrixExpr<ExprType::Mul, M&&, U&&>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] auto operator*(U&& x, M&& m) noexcept requires(!CUDA<M>) {
        return std::forward<M>(m) * std::forward<U>(x);
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard]] auto hadamard(M1&& m1, M2&& m2) noexcept requires(!CUDA<M1> && !CUDA<M2>) {
        return MatrixExpr<ExprType::Mul, M1&&, M2&&>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}

#include "GEVM/Mul.h"
