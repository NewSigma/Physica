/*
 * Copyright 2026 Weibo He.
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
#include "TensorExpr.h"

namespace Physica {
    template<ExprID ID, class LHS, class RHS>
    class device_obj<BinaryTensorExpr<ID, LHS, RHS>> : public device_obj<RValueTensor<TensorExpr<ID, LHS, RHS>>> {
        static_assert(Tensor<LHS> || Tensor<RHS>, "[Error]: Either type should be Tensor");

        using host_obj = BinaryTensorExpr<ID, LHS, RHS>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueTensor<TensorExpr<ID, LHS, RHS>>>;

        using LHS1 = std::remove_reference_t<LHS>;
        using RHS1 = std::remove_reference_t<RHS>;
        using LHS2 = std::conditional<Scalar<LHS>, LHS1, add_device_obj_t<LHS1>>::type;
        using RHS2 = std::conditional<Scalar<RHS>, RHS1, add_device_obj_t<RHS1>>::type;
        using Ref1 = std::conditional<Scalar<LHS>, LHS, add_device_obj_t<LHS>>::type;
        using Ref2 = std::conditional<Scalar<RHS>, RHS, add_device_obj_t<RHS>>::type;
    public:
        using typename Base::IndexType;
    private:
        PlainStruct<LHS2> lhs;
        PlainStruct<RHS2> rhs;
    public:
        __host__ __device__ device_obj(Ref1 lhs_, Ref2 rhs_) noexcept;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ size_t dim(int index) const noexcept;
        [[nodiscard]] __host__ __device__ IndexType getShape() const noexcept;
        [[nodiscard]] __host__ __device__ auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ constexpr static ExprID getExprID() noexcept { return ID; }
    };

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ device_obj<BinaryTensorExpr<ID, LHS, RHS>>::device_obj(Ref1 lhs_, Ref2 rhs_) noexcept
            : lhs(asStruct(lhs_)), rhs(asStruct(rhs_)) {
        if constexpr (Tensor<LHS> && Tensor<RHS>)
            assert(getLHS().getShape() == getRHS().getShape());
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ size_t device_obj<BinaryTensorExpr<ID, LHS, RHS>>::dim(int index) const noexcept {
        if constexpr (Tensor<LHS>)
            return getLHS().dim(index);
        else
            return getRHS().dim(index);
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ auto device_obj<BinaryTensorExpr<ID, LHS, RHS>>::getShape() const noexcept -> IndexType {
        if constexpr (Tensor<LHS>)
            return getLHS().getShape();
        else
            return getRHS().getShape();
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ auto&& device_obj<BinaryTensorExpr<ID, LHS, RHS>>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref1>(self.lhs.getDerived());
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ auto&& device_obj<BinaryTensorExpr<ID, LHS, RHS>>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref2>(self.rhs.getDerived());
    }
}

namespace Physica {
    template<ExprID ID, class Expr1, class Expr2>
    class Traits<device_obj<TensorExpr<ID, Expr1, Expr2>>> : public Traits<TensorExpr<ID, Expr1, Expr2>> {};
}

#include "TensorExprImpl/Operator/Add.cuh"
#include "TensorExprImpl/Operator/Mul.cuh"
