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
        LazyDestroy<M> expr;
    public:
        UnitaryMatrixExpr(M expr_) : expr(std::forward<M>(expr_)) {}
        UnitaryMatrixExpr(const This&) = default;
        UnitaryMatrixExpr(This&&) noexcept requires(isReverseDiff) = default;
        ~UnitaryMatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] decltype(auto) transpose() const noexcept;
        [[nodiscard]] decltype(auto) hermite() const noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return getExpr().getRow(); }
        [[nodiscard]] size_t getCol() const { return getExpr().getCol(); }
        [[nodiscard]] const auto& getExpr() const noexcept { return expr; }
        [[nodiscard]] auto& getExpr() noexcept { return expr; }
    };

    template<ExprType Type, Matrix M>
    decltype(auto) UnitaryMatrixExpr<Type, M>::transpose() const noexcept {
        if constexpr (MatrixOption::isSymmMatrix<M>())
            return Base::getDerived();
        else
            return Base::transpose();
    }

    template<ExprType Type, Matrix M>
    decltype(auto) UnitaryMatrixExpr<Type, M>::hermite() const noexcept {
        if constexpr (MatrixOption::isHermiteMatrix<M>())
            return Base::getDerived();
        else
            return Base::hermite();
    }

    template<ExprType Type, class LHS, class RHS>
    class BinaryMatrixExpr : public RValueMatrix<MatrixExpr<Type, LHS, RHS>> {
        static_assert(Matrix<LHS> || Matrix<RHS>, "[Error]: Either types should be Matrix");

        using Derived = MatrixExpr<Type, LHS, RHS>;
        using This = BinaryMatrixExpr<Type, LHS, RHS>;
        using Base = RValueMatrix<Derived>;
    public:
        using Base::isReverseDiff;
    private:
        LazyDestroy<LHS> lhs;
        LazyDestroy<RHS> rhs;
    public:
        BinaryMatrixExpr(LHS lhs_, RHS rhs_);
        BinaryMatrixExpr(const This&) = default;
        BinaryMatrixExpr(This&&) noexcept requires(isReverseDiff) = default;
        ~BinaryMatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] decltype(auto) transpose() const noexcept;
        [[nodiscard]] decltype(auto) hermite() const noexcept;
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
        [[nodiscard]] auto& getLHS() noexcept { return lhs; }
        [[nodiscard]] auto& getRHS() noexcept { return rhs; }
        /* Static members */
        [[nodiscard]] consteval static bool isStaticSymm() noexcept;
        [[nodiscard]] consteval static bool isStaticHermite() noexcept;
    };

    template<ExprType Type, class LHS, class RHS>
    BinaryMatrixExpr<Type, LHS, RHS>::BinaryMatrixExpr(LHS lhs_, RHS rhs_) : lhs(std::forward<LHS>(lhs_)), rhs(std::forward<RHS>(rhs_)) {
        if constexpr (Matrix<LHS> && Matrix<RHS>) {
            assert(getLHS().getRow() == getRHS().getRow());
            assert(getLHS().getCol() == getRHS().getCol());
        }
        else if constexpr (Vector<RHS>)
            assert(getLHS().getCol() == getRHS().getLength());
    }

    template<ExprType Type, class LHS, class RHS>
    decltype(auto) BinaryMatrixExpr<Type, LHS, RHS>::transpose() const noexcept {
        if constexpr (isStaticSymm())
            return Base::getDerived();
        else
            return Base::transpose();
    }

    template<ExprType Type, class LHS, class RHS>
    decltype(auto) BinaryMatrixExpr<Type, LHS, RHS>::hermite() const noexcept {
        if constexpr (isStaticHermite())
            return Base::getDerived();
        else
            return Base::hermite();
    }

    template<ExprType Type, class LHS, class RHS>
    consteval bool BinaryMatrixExpr<Type, LHS, RHS>::isStaticSymm() noexcept {
        if constexpr (Matrix<LHS> && Matrix<RHS>)
            return MatrixOption::isSymmMatrix<LHS>() && MatrixOption::isSymmMatrix<RHS>();
        else if constexpr (Vector<LHS> || Vector<RHS>)
            return false;
        else if constexpr (Scalar<LHS>)
            return MatrixOption::isSymmMatrix<RHS>();
        else {
            static_assert(Scalar<RHS>, "[Error]: Unexpected type");
            return MatrixOption::isSymmMatrix<LHS>();
        }
    }

    template<ExprType Type, class LHS, class RHS>
    consteval bool BinaryMatrixExpr<Type, LHS, RHS>::isStaticHermite() noexcept {
        if constexpr (Matrix<LHS> && Matrix<RHS>)
            return MatrixOption::isHermiteMatrix<LHS>() && MatrixOption::isHermiteMatrix<RHS>();
        else if constexpr (Vector<LHS> || Vector<RHS>)
            return false;
        else if constexpr (Scalar<LHS>)
            return MatrixOption::isHermiteMatrix<RHS>() && !std::remove_cvref_t<LHS>::isComplex;
        else {
            static_assert(Scalar<RHS>, "[Error]: Unexpected type");
            return MatrixOption::isHermiteMatrix<LHS>() && !std::remove_cvref_t<RHS>::isComplex;
        }
    }
}

namespace Physica {
    template<ExprType Type_, Matrix LHS_, Matrix RHS_>
    class Traits<MatrixExpr<Type_, LHS_, RHS_>> {
    public:
        constexpr static ExprType Type = Type_;
        using LHS = LHS_;
        using RHS = RHS_;
    private:
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
        using ResultType = Internal::BinaryScalarOpRtnTy<typename LHS1::ScalarType, typename RHS1::ScalarType>::Type;

        constexpr static bool SameMajor = MatrixOption::isSameMajor<LHS1, RHS1>();
        constexpr static int Major = SameMajor ? MatrixOption::getMajor<LHS1>()
                                               : int(MatrixOption::AnyMajor);
        constexpr static bool IsReal = Type == ExprType::Abs || Type == ExprType::Square;
    public:
        using ScalarType = std::conditional<IsReal, typename ResultType::RealType, ResultType>::type;
        constexpr static int Option = Major;
        constexpr static size_t RowAtCompile = LHS1::RowAtCompile > RHS1::RowAtCompile ? LHS1::RowAtCompile : RHS1::RowAtCompile;
        constexpr static size_t ColAtCompile = LHS1::ColAtCompile > RHS1::ColAtCompile ? LHS1::ColAtCompile : RHS1::ColAtCompile;
        constexpr static size_t SizeAtCompile = LHS1::SizeAtCompile > RHS1::SizeAtCompile ? LHS1::SizeAtCompile : RHS1::SizeAtCompile;

        constexpr static bool FastAssign = false;
    };

    template<ExprType Type_, Matrix LHS_, Vector RHS_>
    class Traits<MatrixExpr<Type_, LHS_, RHS_>> {
    public:
        constexpr static ExprType Type = Type_;
        using LHS = LHS_;
        using RHS = RHS_;
    private:
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

    template<ExprType Type_, Matrix LHS_, Scalar RHS_>
    class Traits<MatrixExpr<Type_, LHS_, RHS_>> {
    public:
        constexpr static ExprType Type = Type_;
        using LHS = LHS_;
        using RHS = RHS_;
    private:
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename LHS1::ScalarType, RHS1>::Type;
        constexpr static int Option = LHS1::Option;
        constexpr static size_t RowAtCompile = LHS1::RowAtCompile;
        constexpr static size_t ColAtCompile = LHS1::ColAtCompile;
        constexpr static size_t SizeAtCompile = LHS1::SizeAtCompile;

        constexpr static bool FastAssign = Traits<LHS1>::FastAssign;
    };

    template<ExprType Type, Scalar LHS, Matrix RHS>
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
#include "MatrixExprImpl/Ln1p.h"
#include "MatrixExprImpl/Exp.h"
#include "MatrixExprImpl/Sin.h"
#include "MatrixExprImpl/Cos.h"
#include "MatrixExprImpl/Tanh.h"
#include "MatrixExprImpl/Sech.h"
