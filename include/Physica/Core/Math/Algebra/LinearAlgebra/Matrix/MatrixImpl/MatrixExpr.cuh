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

#include "MatrixExpr.h"

namespace Physica {
    template<ExprType Type, Matrix M>
    class device_obj<UnitaryMatrixExpr<Type, M>> : public device_obj<RValueMatrix<MatrixExpr<Type, M>>> {
        using host_obj = UnitaryMatrixExpr<Type, M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<MatrixExpr<Type, M>>>;
        using DeviceMatrix = device_obj<M>;
    private:
        union {
            PlainStruct<const DeviceMatrix> value;
            const DeviceMatrix* ptr;
        } expr;
    public:
        __host__ __device__ inline device_obj(const device_obj<M>& expr_) {
            if constexpr (IsHost())
                expr.value = asStruct(expr_);
            else
                expr.ptr = &expr_;
        }
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getRow() const {
            return getExpr<Owner>().template getRow<Owner>();
        }

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getCol() const {
            return getExpr<Owner>().template getCol<Owner>();
        }

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ const DeviceMatrix& getExpr() const {
            if constexpr (IsHost() || Owner == Side::Host)
                return expr.value.getDerived();
            else
                return *expr.ptr;
        }
    };

    template<ExprType Type, Matrix LHS, class RHS>
    class device_obj<BinaryMatrixExpr<Type, LHS, RHS>> : public device_obj<RValueMatrix<MatrixExpr<Type, LHS, RHS>>> {
        using host_obj = BinaryMatrixExpr<Type, LHS, RHS>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<MatrixExpr<Type, LHS, RHS>>>;
        using DeviceLHS = device_obj<LHS>;
        using DeviceRHS = std::conditional<Scalar<RHS>, typename RHS::ScalarType, device_obj<RHS>>::type;
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
            if constexpr (Matrix<RHS>) {
                assert(lhs_.getRow() == rhs_.getRow());
                assert(lhs_.getCol() == rhs_.getCol());
            }

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
        [[nodiscard]] __host__ __device__ size_t getRow() const {
            return getLHS<Owner>().template getRow<Owner>();
        }

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getCol() const {
            return getLHS<Owner>().template getCol<Owner>();
        }

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
    template <ExprType type, Matrix T1, class T2>
    class Traits<device_obj<MatrixExpr<type, T1, T2>>> : public Traits<MatrixExpr<type, T1, T2>> {};
}

#include "MatrixExprImpl/Add.cuh"
#include "MatrixExprImpl/Mul.cuh"
