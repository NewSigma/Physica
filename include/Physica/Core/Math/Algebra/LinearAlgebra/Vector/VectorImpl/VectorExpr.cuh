/*
 * Copyright 2023-2026 Weibo He.
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

#include "Physica/PlainStruct.h"
#include "VectorExpr.h"

namespace Physica {
    template<ExprID ID, Vector V>
    class device_obj<UnitaryVectorExpr<ID, V>> : public device_obj<RValueVector<VectorExpr<ID, V>>> {
        using host_obj = UnitaryVectorExpr<ID, V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<VectorExpr<ID, V>>>;
        using Ref = add_device_obj<V>::type;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> expr;
    public:
        __host__ __device__ device_obj(Ref expr_) noexcept;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return getExpr().getLength(); }
        [[nodiscard]] __host__ __device__ auto&& getExpr(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static ExprID getExprID() noexcept { return ID; }
        [[nodiscard]] __host__ __device__ constexpr static bool isFastPacket() noexcept;
    };

    template<ExprID ID, Vector V>
    __host__ __device__ device_obj<UnitaryVectorExpr<ID, V>>::device_obj(Ref expr_) noexcept : expr(asStruct(expr_)) {}

    template<ExprID ID, Vector V>
    __host__ __device__ auto&& device_obj<UnitaryVectorExpr<ID, V>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.expr.getDerived());
    }

    template<ExprID ID, Vector V>
    __host__ __device__ constexpr bool device_obj<UnitaryVectorExpr<ID, V>>::isFastPacket() noexcept {
        return host_obj::isFastPacket();
    }

    template<ExprID ID, class LHS, class RHS>
    class device_obj<BinaryVectorExpr<ID, LHS, RHS>> : public device_obj<RValueVector<VectorExpr<ID, LHS, RHS>>> {
        static_assert(Vector<LHS> || Vector<RHS>, "[Error]: Either type should be Vector");

        using host_obj = BinaryVectorExpr<ID, LHS, RHS>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<VectorExpr<ID, LHS, RHS>>>;

        using LHS1 = std::remove_reference_t<LHS>;
        using RHS1 = std::remove_reference_t<RHS>;
        using LHS2 = std::conditional<Scalar<LHS>, LHS1, add_device_obj_t<LHS1>>::type;
        using RHS2 = std::conditional<Scalar<RHS>, RHS1, add_device_obj_t<RHS1>>::type;
        using Ref1 = std::conditional<Scalar<LHS>, LHS, add_device_obj_t<LHS>>::type;
        using Ref2 = std::conditional<Scalar<RHS>, RHS, add_device_obj_t<RHS>>::type;
    private:
        Physica::PlainStruct<LHS2> lhs;
        Physica::PlainStruct<RHS2> rhs;
    public:
        __host__ __device__ device_obj(Ref1 lhs_, Ref2 rhs_) noexcept;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return getLHS().getLength(); }
        [[nodiscard]] __host__ __device__ auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static ExprID getExprID() noexcept { return ID; }
        [[nodiscard]] __host__ __device__ constexpr static bool isFastPacket() noexcept;
    };

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ device_obj<BinaryVectorExpr<ID, LHS, RHS>>::device_obj(Ref1 lhs_, Ref2 rhs_) noexcept : lhs(asStruct(lhs_)), rhs(asStruct(rhs_)) {
        if constexpr (Vector<RHS>)
            assert(getLHS().getLength() == getRHS().getLength());
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ auto&& device_obj<BinaryVectorExpr<ID, LHS, RHS>>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), LHS>(self.lhs.getDerived());
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ auto&& device_obj<BinaryVectorExpr<ID, LHS, RHS>>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), RHS>(self.rhs.getDerived());
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ constexpr bool device_obj<BinaryVectorExpr<ID, LHS, RHS>>::isFastPacket() noexcept {
        return host_obj::isFastPacket();
    }
}

namespace Physica {
    template<ExprID ID, class Expr1, class Expr2>
    class Traits<device_obj<VectorExpr<ID, Expr1, Expr2>>> : public Traits<VectorExpr<ID, Expr1, Expr2>> {};
}

#include "VectorExprImpl/Operator/Add.cuh"
#include "VectorExprImpl/Operator/Sub.cuh"
#include "VectorExprImpl/Operator/Mul.cuh"
#include "VectorExprImpl/Operator/Div.cuh"
#include "VectorExprImpl/Operator/Minus.cuh"
#include "VectorExprImpl/Operator/Reciprocal.cuh"
#include "VectorExprImpl/Operator/Sqrt.cuh"
#include "VectorExprImpl/Operator/Relu.cuh"
#include "VectorExprImpl/Operator/Abs.cuh"
#include "VectorExprImpl/Operator/Unit.cuh"
#include "VectorExprImpl/Operator/Square.cuh"
#include "VectorExprImpl/Operator/Ln.cuh"
#include "VectorExprImpl/Operator/Exp.cuh"
#include "VectorExprImpl/Operator/Sin.cuh"
#include "VectorExprImpl/Operator/Cos.cuh"
#include "VectorExprImpl/Operator/Tan.cuh"
#include "VectorExprImpl/Operator/LnCosh.cuh"
#include "VectorExprImpl/Operator/Softmax.cuh"
