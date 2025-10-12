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
    template<class T, class U>
    class MatrixExpr<ExprType::Sub, T, U>
            : public BinaryMatrixExpr<ExprType::Sub, T, U> {
        static_assert(Scalar<T> || Scalar<U>, "[Error]: Either types should be Scalar");

        using Base = BinaryMatrixExpr<ExprType::Sub, T, U>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc(row, col) - Base::getRHS();
            else
                return Base::getLHS() - Base::getRHS().calc(row, col);
        }

        [[nodiscard]] Tv calc_value(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc_value(row, col) - Base::getRHS().value();
            else
                return Base::getLHS().value() - Base::getRHS().calc_value(row, col);
        }
    };

    template<Matrix M1, Matrix M2>
    class MatrixExpr<ExprType::Sub, M1, M2>
            : public BinaryMatrixExpr<ExprType::Sub, M1, M2> {
        using Base = BinaryMatrixExpr<ExprType::Sub, M1, M2>;
        using This = MatrixExpr<ExprType::Sub, M1, M2>;
        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<M1>() && MatrixOption::isSymmMatrix<M2>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<M1>() && MatrixOption::isHermiteMatrix<M2>();
        using TransposeRtnTy = std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = std::conditional<IsHermite, const This&, Hermite<This>>::type;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;

        void reverse(const auto& grad) const noexcept;

        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
    };

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprType::Sub, M1, M2>::calc(size_t row, size_t col) const -> CoDiff<T> {
        return Base::getLHS().calc(row, col) - Base::getRHS().calc(row, col);
    }

    template<Matrix M1, Matrix M2>
    auto MatrixExpr<ExprType::Sub, M1, M2>::calc_value(size_t row, size_t col) const -> Tv {
        return Base::getLHS().calc_value(row, col) - Base::getRHS().calc_value(row, col);
    }

    template<Matrix M1, Matrix M2>
    void MatrixExpr<ExprType::Sub, M1, M2>::reverse(const auto& grad) const noexcept {
        if constexpr (ReverseDiff<M1>)
            Base::getLHS().reverse(grad);
        if constexpr (ReverseDiff<M2>)
            Base::getRHS().reverse(-grad);
    }

    template<Matrix T, Scalar U>
    [[nodiscard]] auto operator-(T&& m, U&& x) noexcept requires(!CUDA<T>) {
        return MatrixExpr<ExprType::Sub, T&&, U&&>(std::forward<T>(m), std::forward<U>(x));
    }

    template<Matrix T, Scalar U>
    [[nodiscard]] auto operator-(U&& x, T&& m) noexcept requires(!CUDA<T>) {
        return MatrixExpr<ExprType::Sub, U&&, T&&>(std::forward<U>(x), std::forward<T>(m));
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard]] auto operator-(M1&& m1, M2&& m2) noexcept requires(!CUDA<M1> && !CUDA<M2>) {
        return MatrixExpr<ExprType::Sub, M1&&, M2&&>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}
