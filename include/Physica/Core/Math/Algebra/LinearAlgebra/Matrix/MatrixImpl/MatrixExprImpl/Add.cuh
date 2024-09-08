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
    template<class MatrixType, class AnyScalar>
    class device_obj<MatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>>
            : public device_obj<BinaryMatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>> {
        using Base = device_obj<BinaryMatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const {
            return ScalarType(Base::template getLHS<Owner>().calc(row, col)) + ScalarType(Base::template getRHS<Owner>());
        }
    };

    template<class MatrixType1, class MatrixType2>
    class device_obj<MatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>>
            : public device_obj<BinaryMatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>> {
        using Base = device_obj<BinaryMatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const {
            return ScalarType(Base::template getLHS<Owner>().calc(row, col)) + ScalarType(Base::template getRHS<Owner>().calc(row, col));
        }
    };

    template<class MatrixType1, class MatrixType2>
    [[nodiscard]] __host__ __device__ inline auto operator+(
            const device_obj<RValueMatrix<MatrixType1>>& mat1, const device_obj<RValueMatrix<MatrixType2>>& mat2) noexcept {
        return device_obj<MatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>>(mat1.getDerived(), mat2.getDerived());
    }

    template<class MatrixType, class ScalarType>
    [[nodiscard]] __host__ __device__ inline auto operator+(
            const device_obj<RValueMatrix<MatrixType>>& mat, const ScalarBase<ScalarType>& s) noexcept {
        return device_obj<MatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<ScalarType>>>(mat.getDerived(), s.getDerived());
    }
}
