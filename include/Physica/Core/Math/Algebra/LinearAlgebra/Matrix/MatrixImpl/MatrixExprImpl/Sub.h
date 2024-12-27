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

namespace Physica::Core {
    template<Matrix T1, Matrix T2>
    class MatrixExpr<ExprType::Sub, T1, T2>
            : public BinaryMatrixExpr<ExprType::Sub, T1, T2> {
        using Base = BinaryMatrixExpr<ExprType::Sub, T1, T2>;
        using This = MatrixExpr<ExprType::Sub, T1, T2>;
        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<T1>() && MatrixOption::isSymmMatrix<T2>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<T1>() && MatrixOption::isHermiteMatrix<T2>();
        using TransposeRtnTy = std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = std::conditional<IsHermite, const This&, Hermite<This>>::type;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return Base::getLHS().calc(row, col) - Base::getRHS().calc(row, col);
        }

        [[nodiscard]] ValueType calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) - Base::getRHS().calc_value(row, col);
        }

        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
    };

    template<Matrix T, Scalar U>
    class MatrixExpr<ExprType::Sub, T, U>
            : public BinaryMatrixExpr<ExprType::Sub, T, U> {
        using Base = BinaryMatrixExpr<ExprType::Sub, T, U>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return Base::getLHS().calc(row, col) - Base::getRHS();
        }

        [[nodiscard]] ValueType calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) - Base::getRHS().value();
        }
    };

    template<Matrix T, Scalar U>
    [[nodiscard]] inline auto operator-(const T& m, const U& x) noexcept {
        return MatrixExpr<ExprType::Sub, T, U>(m, x);
    }

    template<Matrix T1, Matrix T2>
    [[nodiscard]] inline auto operator-(const T1& m1, const T2& m2) noexcept {
        return MatrixExpr<ExprType::Sub, T1, T2>(m1, m2);
    }
}
