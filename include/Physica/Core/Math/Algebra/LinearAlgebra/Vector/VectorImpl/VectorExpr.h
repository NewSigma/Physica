/*
 * Copyright 2020-2024 Weibo He.
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
#include <Physica/Core/MultiPrecision/Real.h>
#include <Physica/Core/MultiPrecision/ExprType.h>

namespace Physica::Core {
    /**
     * \class VectorExpr implements template expression for vectors, which will reduce temporary objects.
     * 
     * Operations defined as \tparam T1 \tparam Type \tparam T2. e.g. vector + scalar, expression * expression
     */
    template<ExprType Type, class T1, class T2 = T1> class VectorExpr;

    template<ExprType Type, Vector V>
    class UnitaryVectorExpr : public RValueVector<VectorExpr<Type, V>> {
        using This = UnitaryVectorExpr<Type, V>;
        using Base = RValueVector<This>;
    private:
        const V& expr;
    public:
        UnitaryVectorExpr(const RValueVector<V>& expr_) : expr(expr_.getDerived()) {}
        UnitaryVectorExpr(const This&) = delete;
        UnitaryVectorExpr(This&&) noexcept = delete;
        ~UnitaryVectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ size_t getLength() const { return getExpr().getLength(); }
        [[nodiscard]] __host__ __device__ const V& getExpr() const noexcept { return expr; }
    };

    template<ExprType Type, Vector LHS, class RHS>
    class BinaryVectorExpr : public RValueVector<VectorExpr<Type, LHS, RHS>> {
        using This = BinaryVectorExpr<Type, LHS, RHS>;
        using Base = RValueVector<This>;
        using LHS1 = LHS;
        using RHS1 = typename std::conditional<Scalar<RHS>, typename RHS::ScalarType, RHS>::type;
    private:
        const LHS1* lhs;
        const RHS1* rhs;
    public:
        BinaryVectorExpr(const LHS& lhs_, const RHS& rhs_) : lhs(&lhs_) {
            if constexpr (Scalar<RHS>)
                rhs = &rhs_.getDerived();
            else {
                rhs = &rhs_;
                assert(lhs->getLength() == rhs->getLength());
            }
        }
        BinaryVectorExpr(const This&) = delete;
        BinaryVectorExpr(This&&) noexcept = delete;
        ~BinaryVectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ size_t getLength() const { return getLHS().getLength(); }
        [[nodiscard]] __host__ __device__ const LHS1& getLHS() const noexcept { return *lhs; }
        [[nodiscard]] __host__ __device__ const RHS1& getRHS() const noexcept { return *rhs; }
    };
}

namespace Physica {
    template<Core::ExprType Type, class Expr1, class Expr2>
    class Traits<Core::VectorExpr<Type, Expr1, Expr2>> {
        constexpr static size_t Size1 = Traits<Expr1>::SizeAtCompile;
        constexpr static size_t Size2 = Traits<Expr2>::SizeAtCompile;
        constexpr static bool FastAssign1 = Traits<Expr1>::FastAssign;
        constexpr static bool FastAssign2 = Traits<Expr2>::FastAssign;
        constexpr static bool FastPacket1 = Traits<Expr1>::FastPacket;
        constexpr static bool FastPacket2 = Traits<Expr2>::FastPacket;
        constexpr static bool IsAddOrSub = Type == Core::ExprType::Add || Type == Core::ExprType::Sub;

        using ScalarType1 = typename Expr1::ScalarType;
        using RealType = typename ScalarType1::RealType;
        using BinaryScalarType = typename Core::Internal::BinaryScalarOpReturnType<ScalarType1, typename Expr2::ScalarType>::Type;
        static_assert(Size1 == Dynamic || Size2 == Dynamic || (Size1 == Size2), "[Error]: Vector dimentions do not match");
    public:
        using ScalarType = typename std::conditional<Type == Core::ExprType::Abs, RealType, BinaryScalarType>::type;
        constexpr static size_t SizeAtCompile = Size1 > Size2 ? Size1 : Size2;
        constexpr static bool FastAssign = IsAddOrSub && (FastAssign1 || FastAssign2);
        constexpr static bool FastPacket = FastPacket1 && FastPacket2;
    };

    template<Core::ExprType Type, class Expr, Scalar T>
    class Traits<Core::VectorExpr<Type, Expr, T>> {
    public:
        using ScalarType = typename Core::Internal::BinaryScalarOpReturnType<typename Expr::ScalarType, T>::Type;
        constexpr static size_t SizeAtCompile = Traits<Expr>::SizeAtCompile;
        constexpr static bool FastAssign = Traits<Expr>::FastAssign;
        constexpr static bool FastPacket = Traits<Expr>::FastPacket;
    };
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
#include "VectorExprImpl/Exp.h"
#include "VectorExprImpl/Pow.h"
#include "VectorExprImpl/Sin.h"
#include "VectorExprImpl/Cos.h"
#include "VectorExprImpl/Cosh.h"
#include "VectorExprImpl/Tanh.h"
#include "VectorExprImpl/Sech.h"
#include "VectorExprImpl/LnCosh.h"
