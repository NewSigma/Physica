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

#include "Physica/Core/Scalar/ExprID.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/Tensor.h"

namespace Physica {
    template<ExprID ID, Tensor LHS, class RHS = LHS>
    class TensorExpr;

    template<ExprID ID, class LHS, class RHS>
    class BinaryTensorExpr : public RValueTensor<TensorExpr<ID, LHS, RHS>> {
        static_assert(Tensor<LHS> || Tensor<RHS>, "[Error]: Either type should be Tensor");

        using This = BinaryTensorExpr<ID, LHS, RHS>;
        using Base = RValueTensor<TensorExpr<ID, LHS, RHS>>;
    public:
        using Base::isReverseDiff;
    private:
        LazyDestroy<LHS> lhs;
        LazyDestroy<RHS> rhs;
    public:
        BinaryTensorExpr(LHS&& lhs_, RHS&& rhs_);
        BinaryTensorExpr(const This&) = delete;
        BinaryTensorExpr(This&&) noexcept requires(isReverseDiff) = default;
        ~BinaryTensorExpr() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] size_t dim(int index) const noexcept;
        [[nodiscard]] decltype(auto) getShape() const noexcept;
        [[nodiscard]] int getDim() const;
        [[nodiscard]] size_t getSize() const noexcept;
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
    };

    template<ExprID ID, class LHS, class RHS>
    BinaryTensorExpr<ID, LHS, RHS>::BinaryTensorExpr(LHS&& lhs_, RHS&& rhs_) : lhs(std::forward<LHS>(lhs_)), rhs(std::forward<RHS>(rhs_)) {
        if constexpr (Tensor<LHS> && Tensor<RHS>)
            assert(lhs->getShape() == rhs->getShape());
    }

    template<ExprID ID, class LHS, class RHS>
    size_t BinaryTensorExpr<ID, LHS, RHS>::dim(int index) const noexcept {
        if constexpr (Tensor<LHS>)
            return getLHS().getShape(index);
        else
            return getRHS().getShape(index);
    }

    template<ExprID ID, class LHS, class RHS>
    decltype(auto) BinaryTensorExpr<ID, LHS, RHS>::getShape() const noexcept {
        if constexpr (Tensor<LHS>)
            return getLHS().getShape();
        else
            return getRHS().getShape();
    }

    template<ExprID ID, class LHS, class RHS>
    int BinaryTensorExpr<ID, LHS, RHS>::getDim() const {
        if constexpr (Tensor<LHS>)
            return getLHS().getDim();
        else
            return getRHS().getDim();
    }

    template<ExprID ID, class LHS, class RHS>
    size_t BinaryTensorExpr<ID, LHS, RHS>::getSize() const noexcept {
        if constexpr (Tensor<LHS>)
            return getLHS().getSize();
        else
            return getRHS().getSize();
    }

    template<ExprID ID, class LHS, class RHS>
    auto&& BinaryTensorExpr<ID, LHS, RHS>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), LHS>(self.lhs);
    }

    template<ExprID ID, class LHS, class RHS>
    auto&& BinaryTensorExpr<ID, LHS, RHS>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), RHS>(self.rhs);
    }
}

namespace Physica {
    template<ExprID ID, Tensor LHS, Tensor RHS>
    class Traits<TensorExpr<ID, LHS, RHS>> {
        constexpr static int NDim1 = Traits<LHS>::NDim;
        constexpr static int NDim2 = Traits<RHS>::NDim;

        using T1 = std::remove_cvref_t<LHS>::ScalarType;
        using T2 = std::remove_cvref_t<RHS>::ScalarType;
        using T12 = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
        static_assert(NDim1 == Dynamic || NDim2 == Dynamic || (NDim1 == NDim2), "[Error]: Tensor dimentions do not match");
    public:
        using ScalarType = std::conditional<ID == ExprID::Abs, typename T1::RealType, T12>::type;
        constexpr static int NDim = NDim1 > NDim2 ? NDim1 : NDim2;
    };

    template<ExprID ID, Tensor LHS, Scalar RHS>
    class Traits<TensorExpr<ID, LHS, RHS>> {
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename std::remove_cvref_t<LHS>::ScalarType, std::remove_cvref_t<RHS>>::Type;
        constexpr static int NDim = Traits<LHS>::NDim;
    };

    template<ExprID ID, Scalar LHS, Tensor RHS>
    class Traits<TensorExpr<ID, LHS, RHS>> : public Traits<TensorExpr<ID, RHS, LHS>> {};
}

#include "TensorExprImpl/Add.h"
#include "TensorExprImpl/Mul.h"
