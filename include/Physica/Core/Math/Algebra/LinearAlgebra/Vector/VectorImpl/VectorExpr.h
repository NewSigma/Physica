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
#include "Physica/Core/Scalar/ExprID.h"
#include "Physica/Core/Parallel/Parallel.h"
#include "RValueVector.h"

namespace Physica {
    /**
     * \class VectorExpr implements template expression for vectors, which will reduce temporary objects.
     * 
     * Operations defined as \tparam LHS \tparam ID \tparam RHS. e.g. vector + scalar, expression * expression
     */
    template<ExprID ID, class LHS, class RHS = LHS>
    class VectorExpr;

    template<ExprID ID, Vector V>
    class UnitaryVectorExpr : public RValueVector<VectorExpr<ID, V>> {
        static_assert(std::is_reference<V>::value, "[Error]: Expect a reference");
        static_assert(!std::is_const<V>::value, "[Error]: Const is implied");
        using This = UnitaryVectorExpr<ID, V>;
        using Base = RValueVector<VectorExpr<ID, V>>;
    private:
        LazyDestroy<V> expr;
    public:
        UnitaryVectorExpr(V expr_) : expr(std::forward<V>(expr_)) {}
        UnitaryVectorExpr(const This&) = default;
        UnitaryVectorExpr(This&&) noexcept = default;
        ~UnitaryVectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return getExpr().getLength(); }
        [[nodiscard]] const auto& getExpr() const noexcept { return expr; }
        [[nodiscard]] auto& getExpr() noexcept { return expr; }
    };

    template<ExprID ID, class LHS, class RHS>
    class BinaryVectorExpr : public RValueVector<VectorExpr<ID, LHS, RHS>> {
        static_assert(Vector<LHS> || Vector<RHS>, "[Error]: Either type should be Vector");
        static_assert(std::is_reference<LHS>::value, "[Error]: Expect a reference");
        static_assert(std::is_reference<RHS>::value, "[Error]: Expect a reference");
        static_assert(!std::is_const<LHS>::value, "[Error]: Const is implied");
        static_assert(!std::is_const<RHS>::value, "[Error]: Const is implied");

        using This = BinaryVectorExpr<ID, LHS, RHS>;
        using Base = RValueVector<VectorExpr<ID, LHS, RHS>>;
    private:
        LazyDestroy<LHS> lhs;
        LazyDestroy<RHS> rhs;
    public:
        BinaryVectorExpr(LHS lhs_, RHS rhs_);
        BinaryVectorExpr(const This&) = default;
        BinaryVectorExpr(This&&) noexcept = default;
        ~BinaryVectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept;
        [[nodiscard]] const auto& getLHS() const noexcept { return lhs; }
        [[nodiscard]] const auto& getRHS() const noexcept { return rhs; }
        [[nodiscard]] auto& getLHS() noexcept { return lhs; }
        [[nodiscard]] auto& getRHS() noexcept { return rhs; }
    };

    template<ExprID ID, class LHS, class RHS>
    BinaryVectorExpr<ID, LHS, RHS>::BinaryVectorExpr(LHS lhs_, RHS rhs_) : lhs(std::forward<LHS>(lhs_)), rhs(std::forward<RHS>(rhs_)) {
        if constexpr (Vector<LHS> && Vector<RHS>)
            assert(lhs.getLength() == rhs.getLength());
    }

    template<ExprID ID, class LHS, class RHS>
    size_t BinaryVectorExpr<ID, LHS, RHS>::getLength() const noexcept {
        if constexpr (Vector<LHS>)
            return getLHS().getLength();
        else
            return getRHS().getLength();
    }
}

namespace Physica {
    template<ExprID ID_, Vector LHS_, Vector RHS_>
    class Traits<VectorExpr<ID_, LHS_, RHS_>> {
    public:
        constexpr static ExprID ID = ID_;
        using LHS = LHS_;
        using RHS = RHS_;
    private:
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
        using T = LHS1::ScalarType;
        using Tr = T::RealType;
        using T12 = Internal::BinaryScalarOpRtnTy<T, typename RHS1::ScalarType>::Type;

        constexpr static size_t Size1 = Traits<LHS1>::SizeAtCompile;
        constexpr static size_t Size2 = Traits<RHS1>::SizeAtCompile;
        constexpr static bool FastAssign1 = Traits<LHS1>::FastAssign;
        constexpr static bool FastAssign2 = Traits<RHS1>::FastAssign;
        constexpr static bool FastPacket1 = Traits<LHS1>::FastPacket;
        constexpr static bool FastPacket2 = Traits<RHS1>::FastPacket;
        static_assert(Size1 == Dynamic || Size2 == Dynamic || (Size1 == Size2), "[Error]: Vector dimentions do not match");

        constexpr static bool calcFastAssign() {
            if constexpr (ID == ExprID::Minus)
                return Traits<LHS1>::FastAssign;
            else if constexpr (ID == ExprID::Add || ID == ExprID::Sub)
                return FastAssign1 || FastAssign2;
            else
                return false;
        }
    public:
        using ScalarType = std::conditional<ID == ExprID::Abs, Tr, T12>::type;
        constexpr static size_t SizeAtCompile = Size1 > Size2 ? Size1 : Size2;
        constexpr static bool FastAssign = calcFastAssign();
        constexpr static bool FastPacket = FastPacket1 && FastPacket2;
    };

    template<ExprID ID_, Vector LHS_, Scalar RHS_>
    class Traits<VectorExpr<ID_, LHS_, RHS_>> {
    public:
        constexpr static ExprID ID = ID_;
        using LHS = LHS_;
        using RHS = RHS_;
    private:
        using LHS1 = std::remove_cvref_t<LHS>;
        using RHS1 = std::remove_cvref_t<RHS>;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename LHS1::ScalarType, RHS1>::Type;
        constexpr static size_t SizeAtCompile = Traits<LHS1>::SizeAtCompile;
        constexpr static bool FastAssign = Traits<LHS1>::FastAssign;
        constexpr static bool FastPacket = Traits<LHS1>::FastPacket;
    };

    template<ExprID ID, Scalar LHS, Vector RHS>
    class Traits<VectorExpr<ID, LHS, RHS>> : public Traits<VectorExpr<ID, RHS, LHS>> {};
}

#ifdef PHYSICA_MKL
    #include <mkl_vml.h>
#endif
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
#include "VectorExprImpl/Sigmoid.h"
