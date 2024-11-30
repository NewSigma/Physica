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
#include "Physica/Core/MultiPrecision/ExprType.h"

namespace Physica::Core {
    /**
     * \class VectorExpr implements template expression for vectors, which will reduce temporary objects.
     * 
     * Operations defined as \tparam LHS \tparam Type \tparam RHS. e.g. vector + scalar, expression * expression
     */
    template<ExprType Type, Vector LHS, class RHS = LHS> class VectorExpr;

    template<ExprType Type, Vector V>
    class UnitaryVectorExpr : public RValueVector<VectorExpr<Type, V>> {
        using This = UnitaryVectorExpr<Type, V>;
        using Base = RValueVector<This>;
    private:
        const V& expr;
    public:
        UnitaryVectorExpr(const V& expr_) : expr(expr_) {}
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
        using RHS1 = std::conditional<Scalar<RHS>, typename RHS::ScalarType, RHS>::type;
    private:
        const LHS1* lhs;
        const RHS1* rhs;
    public:
        BinaryVectorExpr(const LHS& lhs_, const RHS& rhs_) : lhs(&lhs_), rhs(&rhs_) {
            if constexpr (!Scalar<RHS>)
                assert(lhs->getLength() == rhs->getLength());
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
    template<ExprType Type, Vector LHS, class RHS>
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
        using BinaryScalarType = Core::Internal::BinaryScalarOpRtnTy<ScalarType1, typename RHS::ScalarType>::Type;
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
        using ScalarType = Core::Internal::BinaryScalarOpRtnTy<typename LHS::ScalarType, RHS>::Type;
        constexpr static size_t SizeAtCompile = Traits<LHS>::SizeAtCompile;
        constexpr static bool FastAssign = Traits<LHS>::FastAssign;
        constexpr static bool FastPacket = Traits<LHS>::FastPacket;
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
#include "VectorExprImpl/Tan.h"
#include "VectorExprImpl/Sec.h"
#include "VectorExprImpl/Cosh.h"
#include "VectorExprImpl/Tanh.h"
#include "VectorExprImpl/Sech.h"
#include "VectorExprImpl/LnCosh.h"
