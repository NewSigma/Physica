/*
 * Copyright 2020-2025 Weibo He.
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

#include <cassert>
#include "Physica/Core/Scalar/ExprType.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Parallel/Executor/SeqExecutor.h"

namespace Physica {
    /**
     * \class VectorExpr implements template expression for vectors, which will reduce temporary objects.
     * 
     * Operations defined as \tparam LHS \tparam Type \tparam RHS. e.g. vector + scalar, expression * expression
     */
    template<ExprType Type, class LHS, class RHS = LHS>
    class VectorExpr;

    template<ExprType Type, Vector V>
    class UnitaryVectorExpr : public RValueVector<VectorExpr<Type, V>> {
        using This = UnitaryVectorExpr<Type, V>;
        using Base = RValueVector<VectorExpr<Type, V>>;
    public:
        using Base::isReverseDiff;
    private:
        const V& expr;
    public:
        UnitaryVectorExpr(const V& expr_) : expr(expr_) {}
        UnitaryVectorExpr(const This&) = delete;
        UnitaryVectorExpr(This&&) noexcept requires(isReverseDiff) = default;
        ~UnitaryVectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] size_t getLength() const { return getExpr().getLength(); }
        [[nodiscard]] const V& getExpr() const noexcept { return expr; }
    };

    template<ExprType Type, class LHS, class RHS>
    class BinaryVectorExpr : public RValueVector<VectorExpr<Type, LHS, RHS>> {
        static_assert(Vector<LHS> || Vector<RHS>, "[Error]: Either type should be Vector");

        using This = BinaryVectorExpr<Type, LHS, RHS>;
        using Base = RValueVector<VectorExpr<Type, LHS, RHS>>;
    public:
        using Base::isReverseDiff;
    private:
        const LHS* lhs;
        const RHS* rhs;
    public:
        BinaryVectorExpr(const LHS& lhs_, const RHS& rhs_) : lhs(&lhs_), rhs(&rhs_) {
            if constexpr (Vector<LHS> && Vector<RHS>)
                assert(lhs->getLength() == rhs->getLength());
        }
        BinaryVectorExpr(const This&) = delete;
        BinaryVectorExpr(This&&) noexcept requires(isReverseDiff) = default;
        ~BinaryVectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] size_t getLength() const {
            if constexpr (Vector<LHS>)
                return getLHS().getLength();
            else
                return getRHS().getLength();
        }
        [[nodiscard]] const LHS& getLHS() const noexcept { return *lhs; }
        [[nodiscard]] const RHS& getRHS() const noexcept { return *rhs; }
    };
}

namespace Physica {
    template<ExprType Type, Vector LHS, Vector RHS>
    class Traits<VectorExpr<Type, LHS, RHS>> {
        constexpr static size_t Size1 = Traits<LHS>::SizeAtCompile;
        constexpr static size_t Size2 = Traits<RHS>::SizeAtCompile;
        constexpr static bool FastAssign1 = Traits<LHS>::FastAssign;
        constexpr static bool FastAssign2 = Traits<RHS>::FastAssign;
        constexpr static bool FastPacket1 = Traits<LHS>::FastPacket;
        constexpr static bool FastPacket2 = Traits<RHS>::FastPacket;
        constexpr static bool IsAddOrSub = Type == ExprType::Add || Type == ExprType::Sub;

        using ScalarType1 = LHS::ScalarType;
        using RealType = ScalarType1::RealType;
        using BinaryScalarType = Internal::BinaryScalarOpRtnTy<ScalarType1, typename RHS::ScalarType>::Type;
        static_assert(Size1 == Dynamic || Size2 == Dynamic || (Size1 == Size2), "[Error]: Vector dimentions do not match");
    public:
        using ScalarType = std::conditional<Type == ExprType::Abs, RealType, BinaryScalarType>::type;
        constexpr static size_t SizeAtCompile = Size1 > Size2 ? Size1 : Size2;
        constexpr static bool FastAssign = IsAddOrSub && (FastAssign1 || FastAssign2);
        constexpr static bool FastPacket = FastPacket1 && FastPacket2;
    };

    template<ExprType Type, Vector LHS, Scalar RHS>
    class Traits<VectorExpr<Type, LHS, RHS>> {
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename LHS::ScalarType, RHS>::Type;
        constexpr static size_t SizeAtCompile = Traits<LHS>::SizeAtCompile;
        constexpr static bool FastAssign = Traits<LHS>::FastAssign;
        constexpr static bool FastPacket = Traits<LHS>::FastPacket;
    };

    template<ExprType Type, Scalar LHS, Vector RHS>
    class Traits<VectorExpr<Type, LHS, RHS>> : public Traits<VectorExpr<Type, RHS, LHS>> {};
}

#include "VectorExprImpl/Add.h"
#include "VectorExprImpl/Sub.h"
#include "VectorExprImpl/Mul.h"
#include "VectorExprImpl/Div.h"
#include "VectorExprImpl/Minus.h"
#include "VectorExprImpl/More.h"
#include "VectorExprImpl/MoreEq.h"
#include "VectorExprImpl/Reciprocal.h"
#include "VectorExprImpl/Sqrt.h"
#include "VectorExprImpl/Cbrt.h"
#include "VectorExprImpl/Relu.h"
#include "VectorExprImpl/Abs.h"
#include "VectorExprImpl/Unit.h"
#include "VectorExprImpl/Square.h"
#include "VectorExprImpl/Ln.h"
#include "VectorExprImpl/Ln1p.h"
#include "VectorExprImpl/Ln1pExp.h"
#include "VectorExprImpl/Exp.h"
#include "VectorExprImpl/Pow.h"
#include "VectorExprImpl/Sin.h"
#include "VectorExprImpl/Cos.h"
#include "VectorExprImpl/Tan.h"
#include "VectorExprImpl/Sec.h"
#include "VectorExprImpl/Cosh.h"
#include "VectorExprImpl/Tanh.h"
#include "VectorExprImpl/Sech.h"
#include "VectorExprImpl/ArcSinh.h"
#include "VectorExprImpl/ArcTanh.h"
#include "VectorExprImpl/LnCosh.h"
#include "VectorExprImpl/Softmax.h"
