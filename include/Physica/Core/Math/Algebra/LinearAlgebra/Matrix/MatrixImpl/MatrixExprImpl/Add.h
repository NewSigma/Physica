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
    template<class MatrixType1, class MatrixType2>
    class MatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>
            : public BinaryMatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2> {
        using Base = BinaryMatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>;
        using This = MatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>;
        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<MatrixType1>() && MatrixOption::isSymmMatrix<MatrixType2>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<MatrixType1>() && MatrixOption::isHermiteMatrix<MatrixType2>();
        using TransposeRtnTy = typename std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = typename std::conditional<IsHermite, const This&, Hermite<This>>::type;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(Base::getLHS().calc(row, col)) + ScalarType(Base::getRHS().calc(row, col));
        }

        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
    };

    template<class MatrixType, class AnyScalar>
    class MatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>
            : public BinaryMatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>> {
        using Base = BinaryMatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<AnyScalar>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(Base::getLHS().calc(row, col)) + ScalarType(Base::getRHS());
        }
    };

    template<class MatrixType1, class MatrixType2>
    [[nodiscard]] inline auto operator+(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) noexcept {
        return MatrixExpr<ExpressionType::Add, MatrixType1, MatrixType2>(mat1.getDerived(), mat2.getDerived());
    }

    template<class MatrixType, class ScalarType>
    [[nodiscard]] inline auto operator+(const RValueMatrix<MatrixType>& mat, const ScalarBase<ScalarType>& s) noexcept {
        return MatrixExpr<ExpressionType::Add, MatrixType, ScalarBase<ScalarType>>(mat.getDerived(), s.getDerived());
    }
}
