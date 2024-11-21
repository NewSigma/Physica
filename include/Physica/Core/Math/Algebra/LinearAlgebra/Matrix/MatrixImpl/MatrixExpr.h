/*
 * Copyright 2021-2024 Weibo He.
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

#include "Physica/Core/MultiPrecision/ExprType.h"

namespace Physica::Core {
    /**
     * \class MatrixExpr represents \param T \param type \param U. e.g. matrix + scalar, expression * expression
     */
    template<ExprType Type, Matrix T, class U = T> class MatrixExpr;

    template<ExprType Type, Matrix M>
    class UnitaryMatrixExpr : public RValueMatrix<MatrixExpr<Type, M>> {
        using This = UnitaryMatrixExpr<Type, M>;
        using Base = RValueMatrix<This>;
    private:
        const M& expr;
    public:
        UnitaryMatrixExpr(const M& expr_) : expr(expr_) {}
        UnitaryMatrixExpr(const This&) = delete;
        UnitaryMatrixExpr(This&&) noexcept = delete;
        ~UnitaryMatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return getExpr().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return getExpr().getCol(); }
        [[nodiscard]] __host__ __device__ const M& getExpr() const noexcept { return expr; }
    };

    template<ExprType Type, Matrix LHS, class RHS>
    class BinaryMatrixExpr : public RValueMatrix<MatrixExpr<Type, LHS, RHS>> {
        using This = BinaryMatrixExpr<Type, LHS, RHS>;
        using Base = RValueMatrix<This>;
        using LHS1 = LHS;
        using RHS1 = typename std::conditional<Scalar<RHS>, typename RHS::ScalarType, RHS>::type;
    private:
        const LHS1* lhs;
        const RHS1* rhs;
    public:
        BinaryMatrixExpr(const LHS& lhs_, const RHS& rhs_) : lhs(&lhs_) {
            if constexpr (Scalar<RHS>)
                rhs = &rhs_;
            else {
                rhs = &rhs_;
                assert(lhs->getRow() == rhs->getRow());
                assert(lhs->getCol() == rhs->getCol());
            }
        }
        BinaryMatrixExpr(const This&) = delete;
        BinaryMatrixExpr(This&&) noexcept = delete;
        ~BinaryMatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return getLHS().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return getLHS().getCol(); }
        [[nodiscard]] __host__ __device__ const LHS1& getLHS() const noexcept { return *lhs; }
        [[nodiscard]] __host__ __device__ const RHS1& getRHS() const noexcept { return *rhs; }
    };
}

namespace Physica {
    template<Core::ExprType Type, Matrix T, class U>
    class Traits<Core::MatrixExpr<Type, T, U>> {
        using MatrixOption = Core::MatrixOption;
        constexpr static bool SameMajor = MatrixOption::isSameMajor<T, U>();
        constexpr static int Major = SameMajor ? MatrixOption::getMajor<T>()
                                               : int(MatrixOption::AnyMajor);
        constexpr static bool SameStorage = MatrixOption::isSameStorage<T, U>();
        constexpr static int Storage = SameStorage ? MatrixOption::getStorage<T>()
                                                   : int(MatrixOption::AnyStorage);
        constexpr static bool IsReal = Type == ExprType::Abs || Type == ExprType::Square;
        using ResultType = typename Internal::BinaryScalarOpRtnTy<typename T::ScalarType, typename U::ScalarType>::Type;
    public:
        using ScalarType = typename std::conditional<IsReal, typename ResultType::RealType, ResultType>::type;
        constexpr static int Option = Major | Storage;
        // Optimize: T and U may not have same compiling size, for example, T may be fixed size and U may be dynamic
        constexpr static size_t RowAtCompile = T::RowAtCompile;
        constexpr static size_t ColAtCompile = T::ColAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };

    template<Core::ExprType Type, Matrix T, Scalar U>
    class Traits<Core::MatrixExpr<Type, T, U>> {
    public:
        using ScalarType = typename Core::Internal::BinaryScalarOpRtnTy<typename T::ScalarType, U>::Type;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = T::RowAtCompile;
        constexpr static size_t ColAtCompile = T::ColAtCompile;
        constexpr static size_t SizeAtCompile = T::SizeAtCompile;
    };
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
