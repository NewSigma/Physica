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
    class MatrixExpr<ExpressionType::Mul, MatrixType1, MatrixType2>
            : public RValueMatrix<MatrixExpr<ExpressionType::Mul, MatrixType1, MatrixType2>> {
        using This = MatrixExpr<ExpressionType::Mul, MatrixType1, MatrixType2>;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType1& mat1;
        const MatrixType2& mat2;
    public:
        MatrixExpr(const RValueMatrix<MatrixType1>& mat1_, const RValueMatrix<MatrixType2>& mat2_)
                : mat1(mat1_.getDerived()), mat2(mat2_.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return mat1.calc(row, col) * mat2.calc(row, col);
        }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return mat1.getColumn(); }
    };

    template<class MatrixType, class AnyScalar>
    class MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>
            : public RValueMatrix<MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>> {
        using This = MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<AnyScalar>>;
        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<MatrixType>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<MatrixType>();
        using TransposeRtnTy = typename std::conditional<IsSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = typename std::conditional<IsHermite, const This&, Hermite<This>>::type;
    public:
        using Base = RValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& exp;
        const AnyScalar& scalar;
    public:
        MatrixExpr(const RValueMatrix<MatrixType>& exp_, const ScalarBase<AnyScalar>& base)
                : exp(exp_.getDerived()), scalar(base.getDerived()) {}
        MatrixExpr(const This&) = delete;
        MatrixExpr(This&&) noexcept = delete;
        ~MatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const {
            return ScalarType(exp.calc(row, col)) * ScalarType(scalar);
        }

        [[nodiscard]] ScalarType trace() const { return ScalarType(exp.trace()) * ScalarType(scalar); }
        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return exp.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return exp.getColumn(); }
        [[nodiscard]] const MatrixType& getMatrix() const noexcept { return exp; }
        [[nodiscard]] const AnyScalar& getScalar() const noexcept { return scalar; }
    };

    template<class MatrixType, class ScalarType>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>
    operator*(const ScalarBase<ScalarType>& s, const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>(m, s);
    }

    template<class MatrixType, class ScalarType>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>
    operator*(const RValueMatrix<MatrixType>& m, const ScalarBase<ScalarType>& s) noexcept {
        return MatrixExpr<ExpressionType::Mul, MatrixType, ScalarBase<ScalarType>>(m, s);
    }

    template<class MatrixType1, class MatrixType2>
    [[nodiscard]] inline MatrixExpr<ExpressionType::Mul, MatrixType1, MatrixType2>
    hadamard(const RValueMatrix<MatrixType1>& mat1, const RValueMatrix<MatrixType2>& mat2) noexcept {
        return MatrixExpr<ExpressionType::Mul, MatrixType1, MatrixType2>(mat1, mat2);
    }
}
