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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"

namespace Physica {
    /**
     * \class MatrixExpr represents \param T \param type \param U. e.g. matrix + scalar, expression * expression
     */
    template<ExprType Type, class T, class U = T> class MatrixExpr;

    template<ExprType Type, Matrix M>
    class UnitaryMatrixExpr : public RValueMatrix<MatrixExpr<Type, M>> {
        using Derived = MatrixExpr<Type, M>;
        using This = UnitaryMatrixExpr<Type, M>;
        using Base = RValueMatrix<Derived>;
    public:
        using Base::isReverseDiff;
    private:
        using TransposeRtnTy = std::conditional<MatrixOption::isSymmMatrix<M>(), const Derived&, Transpose<Derived>>::type;
        using HermiteRtnTy = std::conditional<MatrixOption::isHermiteMatrix<M>(), const Derived&, Hermite<Derived>>::type;

        const LazyDestroy<M> expr;
    public:
        UnitaryMatrixExpr(M expr_) : expr(std::forward<M>(expr_)) {}
        UnitaryMatrixExpr(const This&) = default;
        UnitaryMatrixExpr(This&&) noexcept requires(isReverseDiff) = default;
        ~UnitaryMatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return Base::getDerived(); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return Base::getDerived(); }
        /* Getters */
        [[nodiscard]] size_t getRow() const { return getExpr().getRow(); }
        [[nodiscard]] size_t getCol() const { return getExpr().getCol(); }
        [[nodiscard]] const auto& getExpr() const noexcept { return expr; }
    };

    template<ExprType Type, class LHS, class RHS>
    class BinaryMatrixExpr : public RValueMatrix<MatrixExpr<Type, LHS, RHS>> {
        static_assert(Matrix<LHS> || Matrix<RHS>, "[Error]: Either types should be Matrix");

        using Derived = MatrixExpr<Type, LHS, RHS>;
        using This = BinaryMatrixExpr<Type, LHS, RHS>;
        using Base = RValueMatrix<Derived>;
    public:
        using Base::isReverseDiff;
    private:
        using TransposeRtnTy = std::conditional<Traits<Derived>::isSymm, const This&, Transpose<This>>::type;
        using HermiteRtnTy = std::conditional<Traits<Derived>::isHermite, const This&, Hermite<This>>::type;

        const LazyDestroy<LHS> lhs;
        const LazyDestroy<RHS> rhs;
    public:
        BinaryMatrixExpr(LHS lhs_, RHS rhs_) : lhs(std::forward<LHS>(lhs_)), rhs(std::forward<RHS>(rhs_)) {
            if constexpr (Matrix<LHS> && Matrix<RHS>) {
                assert(getLHS().getRow() == getRHS().getRow());
                assert(getLHS().getCol() == getRHS().getCol());
            }
            else if constexpr (Vector<RHS>)
                assert(getLHS().getCol() == getRHS().getLength());
        }
        BinaryMatrixExpr(const This&) = default;
        BinaryMatrixExpr(This&&) noexcept requires(isReverseDiff) = default;
        ~BinaryMatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return Base::getDerived(); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return Base::getDerived(); }
        /* Getters */
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
        [[nodiscard]] const auto& getLHS() const noexcept { return lhs; }
        [[nodiscard]] const auto& getRHS() const noexcept { return rhs; }
    };
}

namespace Physica {
    template<ExprType Type, Matrix LHS, Matrix RHS>
    class Traits<MatrixExpr<Type, LHS, RHS>> {
        using This = MatrixExpr<Type, LHS, RHS>;
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
        using ResultType = Internal::BinaryScalarOpRtnTy<typename LHS1::ScalarType, typename RHS1::ScalarType>::Type;

        constexpr static bool SameMajor = MatrixOption::isSameMajor<LHS1, RHS1>();
        constexpr static int Major = SameMajor ? MatrixOption::getMajor<LHS1>()
                                               : int(MatrixOption::AnyMajor);
        constexpr static bool SameStorage = MatrixOption::isSameStorage<LHS1, RHS1>();
        constexpr static int Storage = SameStorage ? MatrixOption::getStorage<LHS1>()
                                                   : int(MatrixOption::AnyStorage);
        constexpr static bool IsReal = Type == ExprType::Abs || Type == ExprType::Square;
    public:
        using ScalarType = std::conditional<IsReal, typename ResultType::RealType, ResultType>::type;
        constexpr static int Option = Major | Storage;
        constexpr static size_t RowAtCompile = LHS1::RowAtCompile > RHS1::RowAtCompile ? LHS1::RowAtCompile : RHS1::RowAtCompile;
        constexpr static size_t ColAtCompile = LHS1::ColAtCompile > RHS1::ColAtCompile ? LHS1::ColAtCompile : RHS1::ColAtCompile;
        constexpr static size_t SizeAtCompile = LHS1::SizeAtCompile > RHS1::SizeAtCompile ? LHS1::SizeAtCompile : RHS1::SizeAtCompile;

        constexpr static bool FastAssign = false;
        constexpr static bool isSymm = MatrixOption::isSymmMatrix<LHS>() && MatrixOption::isSymmMatrix<RHS>();
        constexpr static bool isHermite = MatrixOption::isHermiteMatrix<LHS>() && MatrixOption::isHermiteMatrix<RHS>();
    };

    template<ExprType Type, Matrix LHS, Scalar RHS>
    class Traits<MatrixExpr<Type, LHS, RHS>> {
        using This = MatrixExpr<Type, LHS, RHS>;
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename LHS1::ScalarType, RHS1>::Type;
        constexpr static int Option = LHS1::Option;
        constexpr static size_t RowAtCompile = LHS1::RowAtCompile;
        constexpr static size_t ColAtCompile = LHS1::ColAtCompile;
        constexpr static size_t SizeAtCompile = LHS1::SizeAtCompile;

        constexpr static bool FastAssign = false;
        constexpr static bool isSymm = MatrixOption::isSymmMatrix<LHS>();
        constexpr static bool isHermite = MatrixOption::isHermiteMatrix<LHS>() && !RHS1::isComplex;
    };

    template<ExprType Type, Scalar LHS, Matrix RHS>
    class Traits<MatrixExpr<Type, LHS, RHS>> : public Traits<MatrixExpr<Type, RHS, LHS>> {};

    template<ExprType Type, Matrix LHS, Vector RHS>
    class Traits<MatrixExpr<Type, LHS, RHS>> {
        using This = MatrixExpr<Type, LHS, RHS>;
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename LHS1::ScalarType, typename RHS1::ScalarType>::Type;
        constexpr static int Option = LHS1::Option;
        constexpr static size_t RowAtCompile = LHS1::RowAtCompile;
        constexpr static size_t ColAtCompile = LHS1::ColAtCompile;
        constexpr static size_t SizeAtCompile = LHS1::SizeAtCompile;

        constexpr static bool FastAssign = false;
        constexpr static bool isSymm = false;
        constexpr static bool isHermite = false;
    };

    template<ExprType Type, Vector LHS, Matrix RHS>
    class Traits<MatrixExpr<Type, LHS, RHS>> : public Traits<MatrixExpr<Type, RHS, LHS>> {};
}

#include "MatrixExprImpl/Minus.h"
#include "MatrixExprImpl/Add.h"
#include "MatrixExprImpl/Sub.h"
#include "MatrixExprImpl/Mul.h"
#include "MatrixExprImpl/Div.h"
#include "MatrixExprImpl/Reciprocal.h"
#include "MatrixExprImpl/Sqrt.h"
#include "MatrixExprImpl/Relu.h"
#include "MatrixExprImpl/Abs.h"
#include "MatrixExprImpl/Unit.h"
#include "MatrixExprImpl/Square.h"
#include "MatrixExprImpl/Ln.h"
#include "MatrixExprImpl/Exp.h"
#include "MatrixExprImpl/Sin.h"
#include "MatrixExprImpl/Cos.h"
#include "MatrixExprImpl/Tanh.h"
#include "MatrixExprImpl/Sech.h"
#include "MatrixExprImpl/ExprVecProduct.h"
