/*
 * Copyright 2020-2026 Weibo He.
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
        using This = UnitaryVectorExpr<ID, V>;
        using Base = RValueVector<VectorExpr<ID, V>>;

        template<std::ranges::view Operand> class View;
    protected:
        using typename Base::T;
    private:
        LazyDestroy<V> expr;
    public:
        UnitaryVectorExpr(V expr_) noexcept : expr(std::forward<V>(expr_)) {}
        UnitaryVectorExpr(const This&) = default;
        UnitaryVectorExpr(This&&) noexcept = default;
        ~UnitaryVectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] constexpr auto view() const noexcept;

        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index, size_t count) const noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return getExpr().getLength(); }
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static ExprID getExprID() noexcept { return ID; }
    };

    template<ExprID ID, Vector V>
    constexpr auto UnitaryVectorExpr<ID, V>::view() const noexcept {
        return View(expr.view());
    }

    template<ExprID ID, Vector V>
    template<int Size>
    auto UnitaryVectorExpr<ID, V>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        return (view().begin() + index).template load<Size>();
    }

    template<ExprID ID, Vector V>
    template<int Size>
    auto UnitaryVectorExpr<ID, V>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size> {
        return (view().begin() + index).template load<Size>(count);
    }

    template<ExprID ID, Vector V>
    auto&& UnitaryVectorExpr<ID, V>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.expr);
    }

    template<ExprID ID, class LHS, class RHS>
    class BinaryVectorExpr : public RValueVector<VectorExpr<ID, LHS, RHS>> {
        static_assert(Vector<LHS> || Vector<RHS>, "[Error]: Either type should be Vector");
        using This = BinaryVectorExpr<ID, LHS, RHS>;
        using Base = RValueVector<VectorExpr<ID, LHS, RHS>>;

        template<class OpLHS, class OpRHS> class View;
    protected:
        using typename Base::T;
    private:
        LazyDestroy<LHS> lhs;
        LazyDestroy<RHS> rhs;
    public:
        BinaryVectorExpr(LHS lhs_, RHS rhs_) noexcept;
        BinaryVectorExpr(const This&) = default;
        BinaryVectorExpr(This&&) noexcept = default;
        ~BinaryVectorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] constexpr auto view() const noexcept;

        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] SIMD<T, Size> packet(size_t index, size_t count) const noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept;
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static ExprID getExprID() noexcept { return ID; }
    };

    template<ExprID ID, class LHS, class RHS>
    BinaryVectorExpr<ID, LHS, RHS>::BinaryVectorExpr(LHS lhs_, RHS rhs_) noexcept : lhs(std::forward<LHS>(lhs_)), rhs(std::forward<RHS>(rhs_)) {
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

    template<ExprID ID, class LHS, class RHS>
    constexpr auto BinaryVectorExpr<ID, LHS, RHS>::view() const noexcept {
        if constexpr (Scalar<LHS>)
            return View(getLHS(), getRHS().view());
        else if constexpr (Scalar<RHS>)
            return View(getLHS().view(), getRHS());
        else
            return View(getLHS().view(), getRHS().view());
    }

    template<ExprID ID, class LHS, class RHS>
    template<int Size>
    auto BinaryVectorExpr<ID, LHS, RHS>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        return (view().begin() + index).template load<Size>();
    }

    template<ExprID ID, class LHS, class RHS>
    template<int Size>
    auto BinaryVectorExpr<ID, LHS, RHS>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size> {
        return (view().begin() + index).template load<Size>(count);
    }

    template<ExprID ID, class LHS, class RHS>
    auto&& BinaryVectorExpr<ID, LHS, RHS>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), LHS>(self.lhs);
    }

    template<ExprID ID, class LHS, class RHS>
    auto&& BinaryVectorExpr<ID, LHS, RHS>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), RHS>(self.rhs);
    }
}

namespace Physica {
    template<ExprID ID, Vector LHS_, Vector RHS_>
    class Traits<VectorExpr<ID, LHS_, RHS_>> {
    public:
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
        constexpr static size_t SizeAtCompile = std::max(Size1, Size2);
        constexpr static bool FastAssign = calcFastAssign();
        constexpr static bool FastPacket = FastPacket1 && FastPacket2;
    };

    template<ExprID ID_, Vector LHS_, Scalar RHS_>
    class Traits<VectorExpr<ID_, LHS_, RHS_>> {
    public:
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

#include "VectorExprImpl/UnitaryVectorExprView.h"
#include "VectorExprImpl/BinaryVectorExprView.h"
#ifdef PHYSICA_MKL
    #include <mkl_vml.h>
#endif
#include "VectorExprImpl/Operator/Add.h"
#include "VectorExprImpl/Operator/Sub.h"
#include "VectorExprImpl/Operator/Mul.h"
#include "VectorExprImpl/Operator/Div.h"
#include "VectorExprImpl/Operator/Minus.h"
#include "VectorExprImpl/Operator/Reciprocal.h"
#include "VectorExprImpl/Operator/Sqrt.h"
#include "VectorExprImpl/Operator/Cbrt.h"
#include "VectorExprImpl/Operator/Relu.h"
#include "VectorExprImpl/Operator/Abs.h"
#include "VectorExprImpl/Operator/Unit.h"
#include "VectorExprImpl/Operator/Square.h"
#include "VectorExprImpl/Operator/Ln.h"
#include "VectorExprImpl/Operator/Ln1p.h"
#include "VectorExprImpl/Operator/Ln1pExp.h"
#include "VectorExprImpl/Operator/Exp.h"
#include "VectorExprImpl/Operator/Pow.h"
#include "VectorExprImpl/Operator/Sin.h"
#include "VectorExprImpl/Operator/Cos.h"
#include "VectorExprImpl/Operator/Tan.h"
#include "VectorExprImpl/Operator/Sec.h"
#include "VectorExprImpl/Operator/Cosh.h"
#include "VectorExprImpl/Operator/Tanh.h"
#include "VectorExprImpl/Operator/Sech.h"
#include "VectorExprImpl/Operator/ArcSinh.h"
#include "VectorExprImpl/Operator/ArcTanh.h"
#include "VectorExprImpl/Operator/LnCosh.h"
#include "VectorExprImpl/Operator/Softmax.h"
#include "VectorExprImpl/Operator/Sigmoid.h"
