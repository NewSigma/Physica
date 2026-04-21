/*
 * Copyright 2021-2026 Weibo He.
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

#include "Physica/Core/Scalar/ExprID.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<ExprID ID, Matrix M>
    class UnitaryMatrixExpr : public RValueMatrix<MatrixExpr<ID, M>> {
        using Derived = MatrixExpr<ID, M>;
        using This = UnitaryMatrixExpr<ID, M>;
        using Base = RValueMatrix<Derived>;
    public:
        using Base::isReverseDiff;
    private:
        LazyDestroy<M> expr;
    public:
        UnitaryMatrixExpr(M expr_) noexcept : expr(std::forward<M>(expr_)) {}
        UnitaryMatrixExpr(const This&) = default;
        UnitaryMatrixExpr(This&&) noexcept requires(isReverseDiff()) = default;
        ~UnitaryMatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] decltype(auto) transpose(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) hermite(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return getExpr().getRow(); }
        [[nodiscard]] size_t getCol() const { return getExpr().getCol(); }
        [[nodiscard]] constexpr auto&& getExpr(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static ExprID getExprID() noexcept { return ID; }
    };

    template<ExprID ID, Matrix M>
    decltype(auto) UnitaryMatrixExpr<ID, M>::transpose(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (self.expr.isStaticSymm())
            return std::forward<Self>(self);
        else {
            using X = Base; // FIXME: clang 22 rejects valid
            [[maybe_unused]] auto x = sizeof(X);
            return std::forward<Self>(self).X::transpose();
        }
    }

    template<ExprID ID, Matrix M>
    decltype(auto) UnitaryMatrixExpr<ID, M>::hermite(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (self.expr.isStaticHermite())
            return std::forward<Self>(self);
        else {
            using X = Base; // FIXME: clang 22 rejects valid
            [[maybe_unused]] auto x = sizeof(X);
            return std::forward<Self>(self).X::hermite();
        }
    }

    template<ExprID ID, Matrix M>
    constexpr auto&& UnitaryMatrixExpr<ID, M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.expr);
    }

    template<ExprID ID, class LHS, class RHS>
    class BinaryMatrixExpr : public RValueMatrix<MatrixExpr<ID, LHS, RHS>> {
        static_assert(Matrix<LHS> || Matrix<RHS>, "[Error]: Either types should be Matrix");

        using Derived = MatrixExpr<ID, LHS, RHS>;
        using This = BinaryMatrixExpr<ID, LHS, RHS>;
        using Base = RValueMatrix<Derived>;
    public:
        using Base::RowAtCompile;
        using Base::isReverseDiff;
    private:
        LazyDestroy<LHS> lhs;
        LazyDestroy<RHS> rhs;
    public:
        BinaryMatrixExpr(LHS lhs_, RHS rhs_) noexcept;
        BinaryMatrixExpr(const This&) = default;
        BinaryMatrixExpr(This&&) noexcept requires(isReverseDiff()) = default;
        ~BinaryMatrixExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;

        [[nodiscard]] auto operator*(this auto&&, Vector auto&& v) noexcept requires(RowAtCompile != 1);
        /* Operations */
        [[nodiscard]] decltype(auto) transpose(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) hermite(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const;
        [[nodiscard]] size_t getCol() const;
        [[nodiscard]] constexpr auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] constexpr auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static ExprID getExprID() noexcept { return ID; }
        [[nodiscard]] consteval static bool isStaticSymm() noexcept;
        [[nodiscard]] consteval static bool isStaticHermite() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept;
    };

    template<ExprID ID, class LHS, class RHS>
    BinaryMatrixExpr<ID, LHS, RHS>::BinaryMatrixExpr(LHS lhs_, RHS rhs_) noexcept : lhs(std::forward<LHS>(lhs_)), rhs(std::forward<RHS>(rhs_)) {
        if constexpr (Matrix<LHS> && Matrix<RHS>) {
            assert(getLHS().getRow() == getRHS().getRow());
            assert(getLHS().getCol() == getRHS().getCol());
        }
        else if constexpr (Vector<RHS>)
            assert(getLHS().getCol() == getRHS().getLength());
    }

    template<ExprID ID, class LHS, class RHS>
    auto BinaryMatrixExpr<ID, LHS, RHS>::operator*(this auto&& self, Vector auto&& v) noexcept requires(RowAtCompile != 1) {
        using V = decltype(v);
        if constexpr (ID == ExprID::Mul && Scalar<RHS>)
            return self.getLHS() * (self.getRHS() * std::forward<V>(v));
        else {
            using Self = decltype(self);
            return GEMV<Self, V>(std::forward<Self>(self), std::forward<V>(v));
        }
    }

    template<ExprID ID, class LHS, class RHS>
    decltype(auto) BinaryMatrixExpr<ID, LHS, RHS>::transpose(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isStaticSymm())
            return std::forward<Self>(self);
        else {
            using X = Base; // FIXME: clang 22 rejects valid
            [[maybe_unused]] auto x = sizeof(X);
            return std::forward<Self>(self).X::transpose();
        }
    }

    template<ExprID ID, class LHS, class RHS>
    decltype(auto) BinaryMatrixExpr<ID, LHS, RHS>::hermite(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isStaticHermite())
            return std::forward<Self>(self);
        else {
            using X = Base; // FIXME: clang 22 rejects valid
            [[maybe_unused]] auto x = sizeof(X);
            return std::forward<Self>(self).X::transpose();
        }
    }

    template<ExprID ID, class LHS, class RHS>
    size_t BinaryMatrixExpr<ID, LHS, RHS>::getRow() const {
        if constexpr (Matrix<LHS>)
            return getLHS().getRow();
        else
            return getRHS().getRow();
    }

    template<ExprID ID, class LHS, class RHS>
    size_t BinaryMatrixExpr<ID, LHS, RHS>::getCol() const {
        if constexpr (Matrix<LHS>)
            return getLHS().getCol();
        else
            return getRHS().getCol();
    }

    template<ExprID ID, class LHS, class RHS>
    constexpr auto&& BinaryMatrixExpr<ID, LHS, RHS>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), LHS>(self.lhs);
    }

    template<ExprID ID, class LHS, class RHS>
    constexpr auto&& BinaryMatrixExpr<ID, LHS, RHS>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), RHS>(self.rhs);
    }

    template<ExprID ID, class LHS, class RHS>
    consteval bool BinaryMatrixExpr<ID, LHS, RHS>::isStaticSymm() noexcept {
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
        if constexpr (Matrix<LHS> && Matrix<RHS>)
            return LHS1::isStaticSymm() && RHS1::isStaticSymm();
        else if constexpr (Vector<LHS> || Vector<RHS>)
            return false;
        else if constexpr (Scalar<LHS>)
            return RHS1::isStaticSymm();
        else {
            static_assert(Scalar<RHS>, "[Error]: Unexpected type");
            return LHS1::isStaticSymm();
        }
    }

    template<ExprID ID, class LHS, class RHS>
    consteval bool BinaryMatrixExpr<ID, LHS, RHS>::isStaticHermite() noexcept {
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
        if constexpr (Matrix<LHS> && Matrix<RHS>)
            return LHS1::isStaticHermite() && RHS1::isStaticHermite();
        else if constexpr (Vector<LHS> || Vector<RHS>)
            return false;
        else if constexpr (Scalar<LHS>)
            return RHS1::isStaticHermite() && !LHS1::isComplex();
        else {
            static_assert(Scalar<RHS>, "[Error]: Unexpected type");
            return LHS1::isStaticHermite() && !RHS1::isComplex();
        }
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ consteval bool BinaryMatrixExpr<ID, LHS, RHS>::isFastAssign() noexcept {
        if constexpr (Scalar<LHS>)
            return std::remove_cvref_t<RHS>::isFastAssign();
        if constexpr (Scalar<RHS>)
            return std::remove_cvref_t<LHS>::isFastAssign();
        return false;
    }
}

namespace Physica {
    template<ExprID ID_, Matrix LHS_, Matrix RHS_>
    class Traits<MatrixExpr<ID_, LHS_, RHS_>> {
    public:
        constexpr static ExprID ID = ID_;
        using LHS = LHS_;
        using RHS = RHS_;
    private:
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
        using ResultType = Internal::BinaryScalarOpRtnTy<typename LHS1::ScalarType, typename RHS1::ScalarType>::Type;

        constexpr static bool SameMajor = MatrixMajor::isSameMajor<LHS1, RHS1>();
        constexpr static bool IsReal = ID == ExprID::Abs || ID == ExprID::Square;
    public:
        using ScalarType = std::conditional<IsReal, typename ResultType::RealType, ResultType>::type;
        constexpr static int Major = SameMajor ? MatrixMajor::getMajor<LHS1>() : MatrixMajor::BothMajor;
        constexpr static size_t RowAtCompile = LHS1::RowAtCompile > RHS1::RowAtCompile ? LHS1::RowAtCompile : RHS1::RowAtCompile;
        constexpr static size_t ColAtCompile = LHS1::ColAtCompile > RHS1::ColAtCompile ? LHS1::ColAtCompile : RHS1::ColAtCompile;
        constexpr static size_t SizeAtCompile = LHS1::SizeAtCompile > RHS1::SizeAtCompile ? LHS1::SizeAtCompile : RHS1::SizeAtCompile;
    };

    template<ExprID ID, Matrix LHS_, Vector RHS_>
    class Traits<MatrixExpr<ID, LHS_, RHS_>> {
    public:
        using LHS = LHS_;
        using RHS = RHS_;
    private:
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename LHS1::ScalarType, typename RHS1::ScalarType>::Type;
        constexpr static int Major = LHS1::Major;
        constexpr static size_t RowAtCompile = LHS1::RowAtCompile;
        constexpr static size_t ColAtCompile = LHS1::ColAtCompile;
        constexpr static size_t SizeAtCompile = LHS1::SizeAtCompile;
    };

    template<ExprID ID, Vector LHS, Matrix RHS>
    class Traits<MatrixExpr<ID, LHS, RHS>> : public Traits<MatrixExpr<ID, RHS, LHS>> {};

    template<ExprID ID, Matrix LHS_, Scalar RHS_>
    class Traits<MatrixExpr<ID, LHS_, RHS_>> {
    public:
        using LHS = LHS_;
        using RHS = RHS_;
    private:
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename LHS1::ScalarType, RHS1>::Type;
        constexpr static int Major = LHS1::Major;
        constexpr static size_t RowAtCompile = LHS1::RowAtCompile;
        constexpr static size_t ColAtCompile = LHS1::ColAtCompile;
        constexpr static size_t SizeAtCompile = LHS1::SizeAtCompile;
    };

    template<ExprID ID, Scalar LHS, Matrix RHS>
    class Traits<MatrixExpr<ID, LHS, RHS>> : public Traits<MatrixExpr<ID, RHS, LHS>> {};
}

#include "MatrixExprImpl/Operator/Minus.h"
#include "MatrixExprImpl/Operator/Add.h"
#include "MatrixExprImpl/Operator/Sub.h"
#include "MatrixExprImpl/Operator/Mul.h"
#include "MatrixExprImpl/Operator/Div.h"
#include "MatrixExprImpl/Operator/Reciprocal.h"
#include "MatrixExprImpl/Operator/Sqrt.h"
#include "MatrixExprImpl/Operator/Relu.h"
#include "MatrixExprImpl/Operator/Abs.h"
#include "MatrixExprImpl/Operator/Unit.h"
#include "MatrixExprImpl/Operator/Square.h"
#include "MatrixExprImpl/Operator/Ln.h"
#include "MatrixExprImpl/Operator/Ln1p.h"
#include "MatrixExprImpl/Operator/Exp.h"
#include "MatrixExprImpl/Operator/Sin.h"
#include "MatrixExprImpl/Operator/Cos.h"
#include "MatrixExprImpl/Operator/Tanh.h"
#include "MatrixExprImpl/Operator/Sech.h"

#include "MatrixExprImpl/Operator/MatrixProduct/GEMV.h"
