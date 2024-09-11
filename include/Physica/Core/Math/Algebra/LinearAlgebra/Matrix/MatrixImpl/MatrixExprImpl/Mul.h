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
    class MatrixExpr<ExprType::Mul, MatrixType1, MatrixType2>
            : public BinaryMatrixExpr<ExprType::Mul, MatrixType1, MatrixType2> {
        using Base = BinaryMatrixExpr<ExprType::Mul, MatrixType1, MatrixType2>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return Base::getLHS().calc(row, col) * Base::getRHS().calc(row, col);
        }
    };

    template<class MatrixType, class AnyScalar>
    class MatrixExpr<ExprType::Mul, MatrixType, ScalarBase<AnyScalar>>
            : public BinaryMatrixExpr<ExprType::Mul, MatrixType, ScalarBase<AnyScalar>> {
        using Base = BinaryMatrixExpr<ExprType::Mul, MatrixType, ScalarBase<AnyScalar>>;
        using This = MatrixExpr<ExprType::Mul, MatrixType, ScalarBase<AnyScalar>>;
        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<MatrixType>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<MatrixType>();
        using TransposeRtnTy = typename std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = typename std::conditional<IsHermite, const This&, Hermite<This>>::type;
    public:
        using typename Base::ScalarType;
        using Base::getLHS;
        using Base::getRHS;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(getLHS().calc(row, col)) * ScalarType(getRHS());
        }

        [[nodiscard]] ScalarType trace() const { return ScalarType(getLHS().trace()) * ScalarType(getRHS()); }
        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
    };

    template<class MatrixType, class ScalarType>
    [[nodiscard]] inline auto operator*(const ScalarBase<ScalarType>& s, const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExpr<ExprType::Mul, MatrixType, ScalarBase<ScalarType>>(m.getDerived(), s.getDerived());
    }

    template<class MatrixType, class ScalarType>
    [[nodiscard]] inline auto operator*(const RValueMatrix<MatrixType>& m, const ScalarBase<ScalarType>& s) noexcept {
        return MatrixExpr<ExprType::Mul, MatrixType, ScalarBase<ScalarType>>(m.getDerived(), s.getDerived());
    }

    template<class MatrixType1, class MatrixType2>
    [[nodiscard]] inline auto hadamard(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) noexcept {
        return MatrixExpr<ExprType::Mul, MatrixType1, MatrixType2>(mat1.getDerived(), mat2.getDerived());
    }
}
