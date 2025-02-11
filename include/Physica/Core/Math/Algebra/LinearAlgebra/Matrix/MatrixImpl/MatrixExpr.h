/*
 * Copyright 2021-2025 Weibo He.
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

#include "Physica/Core/Scalar/ExprType.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"

namespace Physica {
    /**
     * \class MatrixExpr represents \param T \param type \param U. e.g. matrix + scalar, expression * expression
     */
    template<ExprType Type, class T, class U = T> class MatrixExpr;

    template<ExprType Type, Matrix M>
    class UnitaryMatrixExpr : public RValueMatrix<MatrixExpr<Type, M>> {
        using This = UnitaryMatrixExpr<Type, M>;
        using Base = RValueMatrix<MatrixExpr<Type, M>>;
    public:
        using Base::isReverseDiff;
    private:
        const M& expr;
    public:
        UnitaryMatrixExpr(const M& expr_) : expr(expr_) {}
        UnitaryMatrixExpr(const This&) = delete;
        UnitaryMatrixExpr(This&&) noexcept requires(isReverseDiff) = default;
        ~UnitaryMatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] size_t getRow() const { return getExpr().getRow(); }
        [[nodiscard]] size_t getCol() const { return getExpr().getCol(); }
        [[nodiscard]] __host__ __device__ const M& getExpr() const noexcept { return expr; }
    };

    template<ExprType Type, class LHS, class RHS>
    class BinaryMatrixExpr : public RValueMatrix<MatrixExpr<Type, LHS, RHS>> {
        static_assert(Matrix<LHS> || Matrix<RHS>, "[Error]: Either types should be Matrix");

        using This = BinaryMatrixExpr<Type, LHS, RHS>;
        using Base = RValueMatrix<MatrixExpr<Type, LHS, RHS>>;
    public:
        using Base::isReverseDiff;
    private:
        const LHS* lhs;
        const RHS* rhs;
    public:
        BinaryMatrixExpr(const LHS& lhs_, const RHS& rhs_) : lhs(&lhs_), rhs(&rhs_) {
            if constexpr (Matrix<LHS> && Matrix<RHS>) {
                assert(lhs->getRow() == rhs->getRow());
                assert(lhs->getCol() == rhs->getCol());
            }
        }
        BinaryMatrixExpr(const This&) = delete;
        BinaryMatrixExpr(This&&) noexcept requires(isReverseDiff) = default;
        ~BinaryMatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] size_t getRow() const {
            if constexpr (Matrix<LHS>)
                return getLHS().getRow();
            else
                return getRHS().getRow();
        }
        [[nodiscard]] size_t getCol() const {
            if constexpr (Matrix<LHS>)
                return getLHS().getCol();
            else
                return getRHS().getCol();
        }
        [[nodiscard]] __host__ __device__ const LHS& getLHS() const noexcept { return *lhs; }
        [[nodiscard]] __host__ __device__ const RHS& getRHS() const noexcept { return *rhs; }
    };
}

namespace Physica {
    template<ExprType Type, Matrix T, class U>
    class Traits<MatrixExpr<Type, T, U>> {
        constexpr static bool SameMajor = MatrixOption::isSameMajor<T, U>();
        constexpr static int Major = SameMajor ? MatrixOption::getMajor<T>()
                                               : int(MatrixOption::AnyMajor);
        constexpr static bool SameStorage = MatrixOption::isSameStorage<T, U>();
        constexpr static int Storage = SameStorage ? MatrixOption::getStorage<T>()
                                                   : int(MatrixOption::AnyStorage);
        constexpr static bool IsReal = Type == ExprType::Abs || Type == ExprType::Square;
        using ResultType = Internal::BinaryScalarOpRtnTy<typename T::ScalarType, typename U::ScalarType>::Type;
    public:
        using ScalarType = std::conditional<IsReal, typename ResultType::RealType, ResultType>::type;
        constexpr static int Option = Major | Storage;
        // Optimize: T and U may not have same compiling size, for example, T may be fixed size and U may be dynamic
        constexpr static size_t RowAtCompile = T::RowAtCompile;
        constexpr static size_t ColAtCompile = T::ColAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };

    template<ExprType Type, Matrix T, Scalar U>
    class Traits<MatrixExpr<Type, T, U>> {
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T::ScalarType, U>::Type;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = T::RowAtCompile;
        constexpr static size_t ColAtCompile = T::ColAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };

    template<ExprType Type, Scalar T, Matrix U>
    class Traits<MatrixExpr<Type, T, U>> : public Traits<MatrixExpr<Type, U, T>> {};
}

#include "MatrixExprImpl/Minus.h"
#include "MatrixExprImpl/Add.h"
#include "MatrixExprImpl/Sub.h"
#include "MatrixExprImpl/Mul.h"
#include "MatrixExprImpl/Div.h"
#include "MatrixExprImpl/Reciprocal.h"
#include "MatrixExprImpl/Sqrt.h"
#include "MatrixExprImpl/Abs.h"
#include "MatrixExprImpl/Square.h"
#include "MatrixExprImpl/Ln.h"
#include "MatrixExprImpl/Exp.h"
#include "MatrixExprImpl/Sin.h"
#include "MatrixExprImpl/Cos.h"
#include "MatrixExprImpl/ExprVecProduct.h"
