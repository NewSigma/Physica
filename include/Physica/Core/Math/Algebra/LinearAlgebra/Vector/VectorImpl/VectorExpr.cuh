/*
 * Copyright 2023-2024 Weibo He.
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

namespace Physica::Core {
    template<ExprType Type, Vector V>
    class device_obj<UnitaryVectorExpr<Type, V>> : public device_obj<RValueVector<VectorExpr<Type, V>>> {
        using host_obj = UnitaryVectorExpr<Type, V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<VectorExpr<Type, V>>>;
        using DeviceVector = device_obj<V>;
    private:
        union {
            PlainStruct<const DeviceVector> value;
            const DeviceVector* ptr;
        } expr;
    public:
        __host__ __device__ inline device_obj(const device_obj<RValueVector<V>>& expr_) {
            if constexpr (IsHost())
                expr.value = asStruct(expr_.getDerived());
            else
                expr.ptr = &expr_.getDerived();
        }
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const {
            return getExpr<Owner>().template getLength<Owner>();
        }

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ const DeviceVector& getExpr() const {
            if constexpr (IsHost() || Owner == Side::Host)
                return expr.value.getDerived();
            else
                return *expr.ptr;
        }
    };

    template<ExprType Type, Vector LHS, class RHS>
    class device_obj<BinaryVectorExpr<Type, LHS, RHS>>
            : public device_obj<RValueVector<VectorExpr<Type, LHS, RHS>>> {
        using host_obj = BinaryVectorExpr<Type, LHS, RHS>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<VectorExpr<Type, LHS, RHS>>>;
        using DeviceLHS = device_obj<LHS>;
        using DeviceRHS = typename std::conditional<Scalar<RHS>, typename RHS::ScalarType, device_obj<RHS>>::type;
    private:
        union {
            PlainStruct<const DeviceLHS> value;
            const DeviceLHS* ptr;
        } lhs;

        union {
            PlainStruct<const DeviceRHS> value;
            const DeviceRHS* ptr;
        } rhs;
    public:
        __host__ __device__ inline device_obj(const DeviceLHS& lhs_, const DeviceRHS& rhs_) {
            if constexpr (Vector<RHS>)
                assert(lhs_.getLength() == rhs_.getLength());
            if constexpr (IsHost()) {
                lhs.value = asStruct(lhs_);
                rhs.value = asStruct(rhs_);
            }
            else {
                lhs.ptr = &lhs_;
                rhs.ptr = &rhs_;
            }
        }
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const { return getLHS<Owner>().template getLength<Owner>(); }

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ const DeviceLHS& getLHS() const noexcept {
            if constexpr (IsHost() || Owner == Side::Host)
                return lhs.value.getDerived();
            else
                return *lhs.ptr;
        }

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ const DeviceRHS& getRHS() const noexcept {
            if constexpr (IsHost() || Owner == Side::Host)
                return rhs.value.getDerived();
            else
                return *rhs.ptr;
        }
    };
}

namespace Physica {
    template<ExprType Type, class Expr1, class Expr2>
    class Traits<Core::device_obj<Core::VectorExpr<Type, Expr1, Expr2>>> : public Traits<Core::VectorExpr<Type, Expr1, Expr2>> {};
}

#include "VectorExprImpl/Add.cuh"
#include "VectorExprImpl/Sub.cuh"
#include "VectorExprImpl/Mul.cuh"
#include "VectorExprImpl/Div.cuh"
#include "VectorExprImpl/Minus.cuh"
#include "VectorExprImpl/Reciprocal.cuh"
#include "VectorExprImpl/Relu.cuh"
#include "VectorExprImpl/Square.cuh"
#include "VectorExprImpl/Exp.cuh"
