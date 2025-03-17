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

#include "MatrixExpr.h"

namespace Physica {
    template<ExprType Type, Matrix M>
    class device_obj<UnitaryMatrixExpr<Type, M>> : public device_obj<RValueMatrix<MatrixExpr<Type, M>>> {
        static_assert(CUDA<M>, "[Error]: Invalid type");
        using host_obj = UnitaryMatrixExpr<Type, M>;
        using Derived = MatrixExpr<Type, M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<Derived>>;
    private:
        using TransposeRtnTy = std::conditional<MatrixOption::isSymmMatrix<M>(), const device_obj<Derived>&, device_obj<Transpose<Derived>>>::type;
        using HermiteRtnTy = std::conditional<MatrixOption::isHermiteMatrix<M>(), const device_obj<Derived>&, device_obj<Hermite<Derived>>>::type;

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
        [[nodiscard]] __host__ __device__ TransposeRtnTy transpose() const noexcept { return Base::getDerived(); }
        [[nodiscard]] __host__ __device__ HermiteRtnTy hermite() const noexcept { return Base::getDerived(); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return getExpr().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return getExpr().getCol(); }
        [[nodiscard]] __host__ __device__ const auto& getExpr() const { return expr.getDerived(); }
    };

    template<ExprType Type, Matrix M>
    __host__ __device__ device_obj<UnitaryMatrixExpr<Type, M>>::device_obj(M expr_) : expr(asStruct(expr_)) {}

    template<ExprType Type, class LHS, class RHS>
    class device_obj<BinaryMatrixExpr<Type, LHS, RHS>> : public device_obj<RValueMatrix<MatrixExpr<Type, LHS, RHS>>> {
        static_assert(Matrix<LHS> || Matrix<RHS>, "[Error]: Either types should be Matrix");
        static_assert((CUDA<LHS> || Scalar<LHS>) && (CUDA<RHS> || Scalar<RHS>), "[Error]: Invalid type");

        using host_obj = BinaryMatrixExpr<Type, LHS, RHS>;
        using Derived = MatrixExpr<Type, LHS, RHS>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<MatrixExpr<Type, LHS, RHS>>>;
    private:
        using TransposeRtnTy = std::conditional<Traits<Derived>::isSymm, const device_obj<Derived>&, device_obj<Transpose<Derived>>>::type;
        using HermiteRtnTy = std::conditional<Traits<Derived>::isHermite, const device_obj<Derived>&, device_obj<Hermite<Derived>>>::type;

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
        [[nodiscard]] __host__ __device__ TransposeRtnTy transpose() const noexcept { return Base::getDerived(); }
        [[nodiscard]] __host__ __device__ HermiteRtnTy hermite() const noexcept { return Base::getDerived(); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const;
        [[nodiscard]] __host__ __device__ size_t getCol() const;
        [[nodiscard]] __host__ __device__ const auto& getLHS() const noexcept { return lhs.getDerived(); }
        [[nodiscard]] __host__ __device__ const auto& getRHS() const noexcept { return rhs.getDerived(); }
    };

    template<ExprType Type, class LHS, class RHS>
    __host__ __device__ inline device_obj<BinaryMatrixExpr<Type, LHS, RHS>>::device_obj(LHS lhs_, RHS rhs_) : lhs(asStruct(lhs_)), rhs(asStruct(rhs_)) {
        if constexpr (Matrix<LHS> && Matrix<RHS>) {
            assert(getLHS().getRow() == getRHS().getRow());
            assert(getLHS().getCol() == getRHS().getCol());
        }
        else if constexpr (Vector<LHS>)
            assert(getLHS().getLength() == getRHS().getRow());
        else if constexpr (Vector<RHS>)
            assert(getLHS().getRow() == getRHS().getLength());
    }

    template<ExprType Type, class LHS, class RHS>
    __host__ __device__ size_t device_obj<BinaryMatrixExpr<Type, LHS, RHS>>::getRow() const {
        if constexpr (Matrix<LHS>)
            return getLHS().getRow();
        else
            return getRHS().getRow();
    }

    template<ExprType Type, class LHS, class RHS>
    __host__ __device__ size_t device_obj<BinaryMatrixExpr<Type, LHS, RHS>>::getCol() const {
        if constexpr (Matrix<LHS>)
            return getLHS().getCol();
        else
            return getRHS().getCol();
    }
}

namespace Physica {
    template <ExprType type, class T1, class T2>
    class Traits<device_obj<MatrixExpr<type, T1, T2>>> : public Traits<MatrixExpr<type, T1, T2>> {};
}

#include "MatrixExprImpl/Add.cuh"
#include "MatrixExprImpl/Sub.cuh"
#include "MatrixExprImpl/Mul.cuh"
#include "MatrixExprImpl/Div.cuh"
#include "MatrixExprImpl/Sqrt.cuh"
#include "MatrixExprImpl/Relu.cuh"
#include "MatrixExprImpl/Unit.cuh"
#include "MatrixExprImpl/Square.cuh"
#include "MatrixExprImpl/Ln.cuh"
#include "MatrixExprImpl/Tanh.cuh"
#include "MatrixExprImpl/Sech.cuh"
