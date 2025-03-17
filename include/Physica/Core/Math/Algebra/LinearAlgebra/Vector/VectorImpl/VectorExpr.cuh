/*
 * Copyright 2023-2025 Weibo He.
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
    template<ExprType Type, Vector V>
    class device_obj<UnitaryVectorExpr<Type, V>> : public device_obj<RValueVector<VectorExpr<Type, V>>> {
        static_assert(CUDA<V>, "[Error]: Invalid type");
        using host_obj = UnitaryVectorExpr<Type, V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<VectorExpr<Type, V>>>;
    private:
        PlainStruct<const std::remove_cvref_t<V>> expr;
    public:
        __host__ __device__ device_obj(V expr_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const { return getExpr().getLength(); }
        [[nodiscard]] __host__ __device__ const auto& getExpr() const { return expr.getDerived(); }
    };

    template<ExprType Type, Vector V>
    __host__ __device__ device_obj<UnitaryVectorExpr<Type, V>>::device_obj(V expr_) : expr(asStruct(expr_)) {}

    template<ExprType Type, Vector LHS, class RHS>
    class device_obj<BinaryVectorExpr<Type, LHS, RHS>> : public device_obj<RValueVector<VectorExpr<Type, LHS, RHS>>> {
        static_assert(CUDA<LHS> && (CUDA<RHS> || Scalar<RHS>), "[Error]: Invalid type");
        using host_obj = BinaryVectorExpr<Type, LHS, RHS>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<VectorExpr<Type, LHS, RHS>>>;
    private:
        PlainStruct<const std::remove_cvref_t<LHS>> lhs;
        PlainStruct<const std::remove_cvref_t<RHS>> rhs;
    public:
        __host__ __device__ device_obj(LHS lhs_, RHS rhs_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const { return getLHS().getLength(); }
        [[nodiscard]] __host__ __device__ const auto& getLHS() const noexcept { return lhs.getDerived(); }
        [[nodiscard]] __host__ __device__ const auto& getRHS() const noexcept { return rhs.getDerived(); }
    };

    template<ExprType Type, Vector LHS, class RHS>
    __host__ __device__ device_obj<BinaryVectorExpr<Type, LHS, RHS>>::device_obj(LHS lhs_, RHS rhs_) : lhs(asStruct(lhs_)), rhs(asStruct(rhs_)) {
        if constexpr (Vector<RHS>)
            assert(getLHS().getLength() == getRHS().getLength());
    }
}

namespace Physica {
    template<ExprType Type, Vector Expr1, class Expr2>
    class Traits<device_obj<VectorExpr<Type, Expr1, Expr2>>> : public Traits<VectorExpr<Type, Expr1, Expr2>> {};
}

#include "VectorExprImpl/Add.cuh"
#include "VectorExprImpl/Sub.cuh"
#include "VectorExprImpl/Mul.cuh"
#include "VectorExprImpl/Div.cuh"
#include "VectorExprImpl/Minus.cuh"
#include "VectorExprImpl/Reciprocal.cuh"
#include "VectorExprImpl/Sqrt.cuh"
#include "VectorExprImpl/Relu.cuh"
#include "VectorExprImpl/Abs.cuh"
#include "VectorExprImpl/Unit.cuh"
#include "VectorExprImpl/Square.cuh"
#include "VectorExprImpl/Ln.cuh"
#include "VectorExprImpl/Exp.cuh"
#include "VectorExprImpl/Sin.cuh"
#include "VectorExprImpl/Cos.cuh"
#include "VectorExprImpl/Tan.cuh"
#include "VectorExprImpl/LnCosh.cuh"
#include "VectorExprImpl/Softmax.cuh"
