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

#include "Physica/Core/Scalar/ExprType.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/Tensor.h"

namespace Physica {
    template<ExprType Type, Tensor LHS, class RHS = LHS>
    class TensorExpr;

    template<ExprType Type, class LHS, class RHS>
    class BinaryTensorExpr : public RValueTensor<TensorExpr<Type, LHS, RHS>> {
        static_assert(Tensor<LHS> || Tensor<RHS>, "[Error]: Either type should be Tensor");

        using This = BinaryTensorExpr<Type, LHS, RHS>;
        using Base = RValueTensor<TensorExpr<Type, LHS, RHS>>;
    public:
        using Base::isReverseDiff;
    private:
        const LHS* lhs;
        const RHS* rhs;
    public:
        BinaryTensorExpr(const LHS& lhs_, const RHS& rhs_) : lhs(&lhs_), rhs(&rhs_) {
            if constexpr (Tensor<LHS> && Tensor<RHS>)
                assert(lhs->getShape() == rhs->getShape());
        }
        BinaryTensorExpr(const This&) = delete;
        BinaryTensorExpr(This&&) noexcept requires(isReverseDiff) = default;
        ~BinaryTensorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] auto getShape() const noexcept;
        [[nodiscard]] size_t getShape(int dim) const;
        [[nodiscard]] int getDim() const;
        [[nodiscard]] size_t getSize() const noexcept;
        [[nodiscard]] const LHS& getLHS() const noexcept { return *lhs; }
        [[nodiscard]] const RHS& getRHS() const noexcept { return *rhs; }
    };

    template<ExprType Type, class LHS, class RHS>
    auto BinaryTensorExpr<Type, LHS, RHS>::getShape() const noexcept {
        if constexpr (Tensor<LHS>)
            return getLHS().getShape();
        else
            return getRHS().getShape();
    }

    template<ExprType Type, class LHS, class RHS>
    size_t BinaryTensorExpr<Type, LHS, RHS>::getShape(int dim) const {
        if constexpr (Tensor<LHS>)
            return getLHS().getShape(dim);
        else
            return getRHS().getShape(dim);
    }

    template<ExprType Type, class LHS, class RHS>
    int BinaryTensorExpr<Type, LHS, RHS>::getDim() const {
        if constexpr (Tensor<LHS>)
            return getLHS().getDim();
        else
            return getRHS().getDim();
    }

    template<ExprType Type, class LHS, class RHS>
    size_t BinaryTensorExpr<Type, LHS, RHS>::getSize() const noexcept {
        if constexpr (Tensor<LHS>)
            return getLHS().getSize();
        else
            return getRHS().getSize();
    }
}

namespace Physica {
    template<ExprType Type, Tensor LHS, Tensor RHS>
    class Traits<TensorExpr<Type, LHS, RHS>> {
        constexpr static int Dim1 = Traits<LHS>::Dim;
        constexpr static int Dim2 = Traits<RHS>::Dim;

        using T = LHS::ScalarType;
        using Tr = T::RealType;
        using T12 = Internal::BinaryScalarOpRtnTy<T, typename RHS::ScalarType>::Type;
        static_assert(Dim1 == Dynamic || Dim2 == Dynamic || (Dim1 == Dim2), "[Error]: Tensor dimentions do not match");
    public:
        using ScalarType = std::conditional<Type == ExprType::Abs, Tr, T12>::type;
        constexpr static int Dim = Dim1 > Dim2 ? Dim1 : Dim2;
    };

    template<ExprType Type, Tensor LHS, Scalar RHS>
    class Traits<TensorExpr<Type, LHS, RHS>> {
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename LHS::ScalarType, RHS>::Type;
        constexpr static int Dim = Traits<LHS>::Dim;
    };

    template<ExprType Type, Scalar LHS, Tensor RHS>
    class Traits<TensorExpr<Type, LHS, RHS>> : public Traits<TensorExpr<Type, RHS, LHS>> {};
}

#include "TensorExprImpl/Add.h"
#include "TensorExprImpl/Mul.h"
