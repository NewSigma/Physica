/*
 * Copyright 2024 Weibo He.
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
    class MatrixExpr<ExprType::Add, M, U>
            : public BinaryMatrixExpr<ExprType::Add, M, U> {
        using Base = BinaryMatrixExpr<ExprType::Add, M, U>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const {
            return Base::getLHS().calc(row, col) + Base::getRHS();
        }

        [[nodiscard]] Tv calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) + Base::getRHS().value();
        }
    };

    template<Matrix M, Vector U>
    class MatrixExpr<ExprType::Add, M, U>
            : public BinaryMatrixExpr<ExprType::Add, M, U> {
        using Base = BinaryMatrixExpr<ExprType::Add, M, U>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const {
            return Base::getLHS().calc(row, col) + Base::getRHS().calc(row);
        }

        [[nodiscard]] Tv calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) + Base::getRHS().calc_value(row);
        }
    };

    template<Matrix T1, Matrix T2>
    class MatrixExpr<ExprType::Add, T1, T2>
            : public BinaryMatrixExpr<ExprType::Add, T1, T2> {
        using Base = BinaryMatrixExpr<ExprType::Add, T1, T2>;
        using This = MatrixExpr<ExprType::Add, T1, T2>;
        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<T1>() && MatrixOption::isSymmMatrix<T2>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<T1>() && MatrixOption::isHermiteMatrix<T2>();
        using TransposeRtnTy = std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = std::conditional<IsHermite, const This&, Hermite<This>>::type;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const {
            return Base::getLHS().calc(row, col) + Base::getRHS().calc(row, col);
        }

        [[nodiscard]] Tv calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) + Base::getRHS().calc_value(row, col);
        }

        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
    };

    template<Matrix M, Scalar U>
    [[nodiscard]] auto operator+(M&& m, U&& x) noexcept requires(!CUDA<M>) {
        return MatrixExpr<ExprType::Add, M&&, U&&>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] auto operator+(U&& x, M&& m) noexcept requires(!CUDA<M>) {
        return m + x;
    }

    template<Matrix M, Vector U>
    [[nodiscard]] auto operator+(M&& m, U&& x) noexcept requires(!CUDA<M> && !CUDA<U>) {
        return MatrixExpr<ExprType::Add, M&&, U&&>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Vector U>
    [[nodiscard]] auto operator+(U&& x, M&& m) noexcept requires(!CUDA<M> && !CUDA<U>) {
        return m + x;
    }

    template<Matrix T1, Matrix T2>
    [[nodiscard]] auto operator+(T1&& m1, T2&& m2) noexcept requires(!CUDA<T1> && !CUDA<T2>) {
        return MatrixExpr<ExprType::Add, T1&&, T2&&>(std::forward<T1>(m1), std::forward<T2>(m2));
    }
}
