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
    template<ExprID ID, Vector V>
    class UnaryVectorExpr : public RValueVector<VectorExpr<ID, V>> {
        using This = UnaryVectorExpr<ID, V>;
        using Base = RValueVector<VectorExpr<ID, V>>;

        template<std::ranges::view Operand> class View;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<V> expr;
    public:
        UnaryVectorExpr(V expr_) noexcept : expr(std::forward<V>(expr_)) {}
        UnaryVectorExpr(const This&) = default;
        UnaryVectorExpr(This&&) noexcept = default;
        ~UnaryVectorExpr() = default;
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
        [[nodiscard]] constexpr auto&& getExpr(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static ExprID getExprID() noexcept { return ID; }
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept;
        [[nodiscard]] __host__ __device__ constexpr static bool isFastPacket() noexcept;
    };

    template<ExprID ID, Vector V>
    constexpr auto UnaryVectorExpr<ID, V>::view() const noexcept {
        return View(expr.view());
    }

    template<ExprID ID, Vector V>
    template<int Size>
    auto UnaryVectorExpr<ID, V>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        return (view().begin() + index).template load<Size>();
    }

    template<ExprID ID, Vector V>
    template<int Size>
    auto UnaryVectorExpr<ID, V>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size> {
        return (view().begin() + index).template load<Size>(count);
    }

    template<ExprID ID, Vector V>
    constexpr auto&& UnaryVectorExpr<ID, V>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.expr);
    }

    template<ExprID ID, Vector V>
    __host__ __device__ consteval size_t UnaryVectorExpr<ID, V>::getSizeAtCompile() noexcept {
        return std::remove_cvref_t<V>::getSizeAtCompile();
    }

    template<ExprID ID, Vector V>
    __host__ __device__ consteval bool UnaryVectorExpr<ID, V>::isFastAssign() noexcept {
        if constexpr (ID == ExprID::Minus)
            return std::remove_cvref_t<V>::isFastAssign();
        return false;
    }

    template<ExprID ID, Vector V>
    __host__ __device__ constexpr bool UnaryVectorExpr<ID, V>::isFastPacket() noexcept {
        return std::remove_cvref_t<V>::isFastPacket();
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
        decay_rvalue_t<LHS> lhs;
        decay_rvalue_t<RHS> rhs;
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
        [[nodiscard]] constexpr auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] constexpr auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static ExprID getExprID() noexcept { return ID; }
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept;
        [[nodiscard]] __host__ __device__ constexpr static bool isFastPacket() noexcept;
    };

    template<ExprID ID, class LHS, class RHS>
    BinaryVectorExpr<ID, LHS, RHS>::BinaryVectorExpr(LHS lhs_, RHS rhs_) noexcept : lhs(std::forward<LHS>(lhs_)), rhs(std::forward<RHS>(rhs_)) {
        if constexpr (Vector<LHS> && Vector<RHS>) {
            constexpr bool Size1 = lhs.getSizeAtCompile();
            constexpr bool Size2 = rhs.getSizeAtCompile();
            static_assert(Size1 == Size2 || Size1 == Dynamic || Size2 == Dynamic);
            assert(lhs.getLength() == rhs.getLength());
        }
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
    constexpr auto&& BinaryVectorExpr<ID, LHS, RHS>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), LHS>(self.lhs);
    }

    template<ExprID ID, class LHS, class RHS>
    constexpr auto&& BinaryVectorExpr<ID, LHS, RHS>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), RHS>(self.rhs);
    }

    template<ExprID ID, class LHS, class RHS>
    consteval size_t BinaryVectorExpr<ID, LHS, RHS>::getSizeAtCompile() noexcept {
        if constexpr (Scalar<LHS>)
            return std::remove_cvref_t<RHS>::getSizeAtCompile();
        else if constexpr (Scalar<RHS>)
            return std::remove_cvref_t<LHS>::getSizeAtCompile();
        else
            return std::max(std::remove_cvref_t<LHS>::getSizeAtCompile(), std::remove_cvref_t<RHS>::getSizeAtCompile());
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ consteval bool BinaryVectorExpr<ID, LHS, RHS>::isFastAssign() noexcept {
        if constexpr (Scalar<LHS>)
            return std::remove_cvref_t<RHS>::isFastAssign();
        else if constexpr (Scalar<RHS>)
            return std::remove_cvref_t<LHS>::isFastAssign();

        if constexpr (Vector<LHS> && Vector<RHS>) {
            constexpr bool FastAssign1 = std::remove_cvref_t<RHS>::isFastAssign();
            constexpr bool FastAssign2 = std::remove_cvref_t<LHS>::isFastAssign();
            if constexpr (ID == ExprID::Add || ID == ExprID::Sub)
                return FastAssign1 || FastAssign2;
        }
        return false;
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ constexpr bool BinaryVectorExpr<ID, LHS, RHS>::isFastPacket() noexcept {
        if constexpr (Scalar<LHS>)
            return std::remove_cvref_t<RHS>::isFastPacket();
        else if constexpr (Scalar<RHS>)
            return std::remove_cvref_t<LHS>::isFastPacket();
        else
            return std::remove_cvref_t<LHS>::isFastPacket() && std::remove_cvref_t<RHS>::isFastPacket();
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
    public:
        using ScalarType = std::conditional<ID == ExprID::Abs, Tr, T12>::type;
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
    };

    template<ExprID ID, Scalar LHS, Vector RHS>
    class Traits<VectorExpr<ID, LHS, RHS>> : public Traits<VectorExpr<ID, RHS, LHS>> {};
}

#include "VectorExprImpl/UnaryVectorExprView.h"
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
#include "VectorExprImpl/Operator/Softplus.h"
