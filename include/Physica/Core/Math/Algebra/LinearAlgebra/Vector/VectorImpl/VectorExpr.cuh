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
    template<ExpressionType ExprType, class VectorType>
    class device_obj<UnitaryVectorExpr<ExprType, VectorType>> : public device_obj<RValueVector<VectorExpr<ExprType, VectorType>>> {
        using host_obj = UnitaryVectorExpr<ExprType, VectorType>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<VectorExpr<ExprType, VectorType>>>;
        using DeviceVector = device_obj<VectorType>;
    private:
        union {
            PlainStruct<const DeviceVector> value;
            const DeviceVector* ptr;
        } expr;
    public:
        __host__ __device__ inline device_obj(const device_obj<RValueVector<VectorType>>& expr_) {
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

    template<ExpressionType ExprType, class LHS, class RHS>
    class device_obj<BinaryVectorExpr<ExprType, LHS, RHS>>
            : public device_obj<RValueVector<VectorExpr<ExprType, LHS, RHS>>> {
        static_assert(Internal::is_vector<LHS>::value, "[Error]: Invalid left hand side type");
        using host_obj = BinaryVectorExpr<ExprType, LHS, RHS>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<VectorExpr<ExprType, LHS, RHS>>>;
        using DeviceLHS = device_obj<LHS>;
        using DeviceRHS = typename std::conditional<is_scalar<RHS>::value, typename RHS::ScalarType, device_obj<RHS>>::type;
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
            if constexpr (Internal::is_vector<RHS>::value)
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
    template<ExpressionType Type, class Exp1, class Exp2>
    class Traits<Core::device_obj<Core::VectorExpr<Type, Exp1, Exp2>>> : public Traits<Core::VectorExpr<Type, Exp1, Exp2>> {};
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
