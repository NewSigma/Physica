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

#include "MatrixExpr.h"

namespace Physica {
    template<ExprID ID, Matrix M>
    class device_obj<UnitaryMatrixExpr<ID, M>> : public device_obj<RValueMatrix<MatrixExpr<ID, M>>> {
        static_assert(CUDA<M>, "[Error]: Invalid type");
        using host_obj = UnitaryMatrixExpr<ID, M>;
        using Derived = MatrixExpr<ID, M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<Derived>>;
    private:
        PlainStruct<const std::remove_cvref_t<M>> expr;
    public:
        __host__ __device__ device_obj(M expr_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ decltype(auto) transpose() const noexcept;
        [[nodiscard]] __host__ __device__ decltype(auto) hermite() const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return getExpr().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return getExpr().getCol(); }
        [[nodiscard]] __host__ __device__ auto&& getExpr(this auto&&) noexcept;
    };

    template<ExprID ID, Matrix M>
    __host__ __device__ device_obj<UnitaryMatrixExpr<ID, M>>::device_obj(M expr_) : expr(asStruct(expr_)) {}

    template<ExprID ID, Matrix M>
    __host__ __device__ decltype(auto) device_obj<UnitaryMatrixExpr<ID, M>>::transpose() const noexcept {
        if constexpr (MatrixOption::isSymmMatrix<M>())
            return Base::getDerived();
        else
            return Base::transpose();
    }

    template<ExprID ID, Matrix M>
    __host__ __device__ decltype(auto) device_obj<UnitaryMatrixExpr<ID, M>>::hermite() const noexcept {
        if constexpr (MatrixOption::isHermiteMatrix<M>())
            return Base::getDerived();
        else
            return Base::hermite();
    }

    template<ExprID ID, Matrix M>
    __host__ __device__ auto&& device_obj<UnitaryMatrixExpr<ID, M>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.expr.getDerived());
    }

    template<ExprID ID, class LHS, class RHS>
    class device_obj<BinaryMatrixExpr<ID, LHS, RHS>> : public device_obj<RValueMatrix<MatrixExpr<ID, LHS, RHS>>> {
        static_assert(Matrix<LHS> || Matrix<RHS>, "[Error]: Either types should be Matrix");
        static_assert((CUDA<LHS> || Scalar<LHS>) && (CUDA<RHS> || Scalar<RHS>), "[Error]: Invalid type");

        using host_obj = BinaryMatrixExpr<ID, LHS, RHS>;
        using Derived = MatrixExpr<ID, LHS, RHS>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<MatrixExpr<ID, LHS, RHS>>>;
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
        /* Operations */
        [[nodiscard]] __host__ __device__ decltype(auto) transpose() const noexcept;
        [[nodiscard]] __host__ __device__ decltype(auto) hermite() const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const;
        [[nodiscard]] __host__ __device__ size_t getCol() const;
        [[nodiscard]] __host__ __device__ auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto&& getRHS(this auto&&) noexcept;
    };

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ device_obj<BinaryMatrixExpr<ID, LHS, RHS>>::device_obj(LHS lhs_, RHS rhs_) : lhs(asStruct(lhs_)), rhs(asStruct(rhs_)) {
        if constexpr (Matrix<LHS> && Matrix<RHS>) {
            assert(getLHS().getRow() == getRHS().getRow());
            assert(getLHS().getCol() == getRHS().getCol());
        }
        else if constexpr (Vector<LHS>)
            assert(getLHS().getLength() == getRHS().getRow());
        else if constexpr (Vector<RHS>)
            assert(getLHS().getRow() == getRHS().getLength());
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ decltype(auto) device_obj<BinaryMatrixExpr<ID, LHS, RHS>>::transpose() const noexcept {
        if constexpr (host_obj::isStaticSymm())
            return Base::getDerived();
        else
            return Base::transpose();
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ decltype(auto) device_obj<BinaryMatrixExpr<ID, LHS, RHS>>::hermite() const noexcept {
        if constexpr (host_obj::isStaticHermite())
            return Base::getDerived();
        else
            return Base::hermite();
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ size_t device_obj<BinaryMatrixExpr<ID, LHS, RHS>>::getRow() const {
        if constexpr (Matrix<LHS>)
            return getLHS().getRow();
        else
            return getRHS().getRow();
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ size_t device_obj<BinaryMatrixExpr<ID, LHS, RHS>>::getCol() const {
        if constexpr (Matrix<LHS>)
            return getLHS().getCol();
        else
            return getRHS().getCol();
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ auto&& device_obj<BinaryMatrixExpr<ID, LHS, RHS>>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), LHS>(self.lhs.getDerived());
    }

    template<ExprID ID, class LHS, class RHS>
    __host__ __device__ auto&& device_obj<BinaryMatrixExpr<ID, LHS, RHS>>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), RHS>(self.rhs.getDerived());
    }
}

namespace Physica {
    template <ExprID ID, class T1, class T2>
    class Traits<device_obj<MatrixExpr<ID, T1, T2>>> : public Traits<MatrixExpr<ID, T1, T2>> {};
}

#include "MatrixExprImpl/Minus.cuh"
#include "MatrixExprImpl/Add.cuh"
#include "MatrixExprImpl/Sub.cuh"
#include "MatrixExprImpl/Mul.cuh"
#include "MatrixExprImpl/Div.cuh"
#include "MatrixExprImpl/Sqrt.cuh"
#include "MatrixExprImpl/Relu.cuh"
#include "MatrixExprImpl/Unit.cuh"
#include "MatrixExprImpl/Square.cuh"
#include "MatrixExprImpl/Ln.cuh"
#include "MatrixExprImpl/Sin.cuh"
#include "MatrixExprImpl/Cos.cuh"
#include "MatrixExprImpl/Cosh.cuh"
#include "MatrixExprImpl/Sinh.cuh"
#include "MatrixExprImpl/Tanh.cuh"
#include "MatrixExprImpl/Sech.cuh"
