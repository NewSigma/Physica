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
    template<class MatrixType, Scalar T>
    class device_obj<MatrixExpr<ExprType::Mul, MatrixType, T>>
            : public device_obj<BinaryMatrixExpr<ExprType::Mul, MatrixType, T>> {
        using Base = device_obj<BinaryMatrixExpr<ExprType::Mul, MatrixType, T>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const {
            return ScalarType(Base::template getLHS<Owner>().calc(row, col)) * ScalarType(Base::template getRHS<Owner>());
        }
    };

    template<class MatrixType, Scalar T>
    [[nodiscard]] __device__ inline auto operator*(const T& x, const device_obj<RValueMatrix<MatrixType>>& m) noexcept {
        return device_obj<MatrixExpr<ExprType::Mul, MatrixType, T>>(m.getDerived(), x);
    }

    template<class MatrixType, Scalar T>
    [[nodiscard]] __device__ inline auto operator*(const device_obj<RValueMatrix<MatrixType>>& m, const T& x) noexcept {
        return device_obj<MatrixExpr<ExprType::Mul, MatrixType, T>>(m.getDerived(), x);
    }
}
