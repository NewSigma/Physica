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
        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<M>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<M>();
        using TransposeRtnTy = std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = std::conditional<IsHermite, const This&, Hermite<This>>::type;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operators */
        [[nodiscard]] auto operator-() const& noexcept;
        [[nodiscard]] auto operator-() && noexcept;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        [[nodiscard]] T trace() const { return getLHS().trace() * getRHS(); }
        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
        [[nodiscard]] T sum() const { return Base::getLHS().sum() * Base::getRHS(); }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprType::Mul, M, U>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return getLHS().calc(row, col) * getRHS();
    }

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprType::Mul, M, U>::calc_value(size_t row, size_t col) const -> Tv {
        return getLHS().calc_value(row, col) * getRHS().value();
    }

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprType::Mul, M, U>::operator-() const& noexcept {
        return getLHS() * (-getRHS());
    }

    template<Matrix M, Scalar U>
    auto MatrixExpr<ExprType::Mul, M, U>::operator-() && noexcept {
        return std::move(getLHS()) * (-getRHS());
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
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const {
            return getLHS().calc(row, col) * getRHS().calc(row, col);
        }

        [[nodiscard]] Tv calc_value(size_t row, size_t col) const {
            return getLHS().calc_value(row, col) * getRHS().calc_value(row, col);
        }

        void reverse(const Matrix auto& grad) const noexcept requires(isReverseDiff);
        using Base::reverse;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprType::Mul, M1, M2>::reverse(const Matrix auto& grad) const noexcept requires(isReverseDiff) {
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        if constexpr (Diffable<M1>)
            lhs.reverse(hadamard(rhs.values(), grad));
        if constexpr (Diffable<M2>)
            rhs.reverse(hadamard(lhs.values(), grad));
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] inline auto operator*(M&& m, U&& x) noexcept requires(!CUDA<M>) {
        return MatrixExpr<ExprType::Mul, M&&, U&&>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] inline auto operator*(U&& x, M&& m) noexcept requires(!CUDA<M>) {
        return m * x;
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard]] inline auto hadamard(M1&& m1, M2&& m2) noexcept requires(!CUDA<M1> && !CUDA<M2>) {
        return MatrixExpr<ExprType::Mul, M1&&, M2&&>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}
