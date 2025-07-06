/*
 * Copyright 2024-2025 Weibo He.
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

#include "../MatrixExpr.cuh"

namespace Physica {
    template<Matrix M, Scalar U>
    class device_obj<MatrixExpr<ExprType::Mul, M, U>>
            : public device_obj<BinaryMatrixExpr<ExprType::Mul, M, U>> {
        using Base = device_obj<BinaryMatrixExpr<ExprType::Mul, M, U>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const {
            if constexpr (isReverseDiff)
                return calc_value(row, col);
            else
                return Base::getLHS().calc(row, col) * Base::getRHS();
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) * Base::getRHS().value();
        }
    };

    template<Matrix T1, Vector T2>
    class device_obj<MatrixExpr<ExprType::Mul, T1, T2>>
            : public device_obj<BinaryMatrixExpr<ExprType::Mul, T1, T2>> {
        using Base = device_obj<BinaryMatrixExpr<ExprType::Mul, T1, T2>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const {
            if constexpr (isReverseDiff)
                return calc_value(row, col);
            else
                return Base::getLHS().calc(row, col) * Base::getRHS().calc(row);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) * Base::getRHS().calc_value(row);
        }
    };

    template<Matrix M1, Matrix M2>
    class device_obj<MatrixExpr<ExprType::Mul, M1, M2>>
            : public device_obj<BinaryMatrixExpr<ExprType::Mul, M1, M2>> {
        using Base = device_obj<BinaryMatrixExpr<ExprType::Mul, M1, M2>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const {
            if constexpr (isReverseDiff)
                return calc_value(row, col);
            else
                return getLHS().calc(row, col) * getRHS().calc(row, col);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            return getLHS().calc_value(row, col) * getRHS().calc_value(row, col);
        }

        void reverse(const Matrix auto& grad) const noexcept requires(isReverseDiff);
        using Base::reverse;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Matrix M1, Matrix M2>
    void device_obj<MatrixExpr<ExprType::Mul, M1, M2>>::reverse(const Matrix auto& grad) const noexcept requires(isReverseDiff) {
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        if constexpr (Diffable<M1>)
            lhs.reverse(hadamard(rhs.values(), grad));
        if constexpr (Diffable<M2>)
            rhs.reverse(hadamard(lhs.values(), grad));
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] __host__ __device__ inline auto operator*(M&& m, U&& x) noexcept requires(CUDA<M>) {
        return device_obj<MatrixExpr<ExprType::Mul, M&&, U&&>>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] __host__ __device__ inline auto operator*(U&& x, M&& m) noexcept requires(CUDA<M>) {
        return m * x;
    }

    template<Matrix M, Vector V>
    [[nodiscard]] __host__ __device__ inline auto hadamard(M&& m, V&& x) noexcept requires(CUDA<M> && CUDA<V>) {
        return device_obj<MatrixExpr<ExprType::Mul, M&&, V&&>>(std::forward<M>(m), std::forward<V>(x));
    }

    template<Matrix M, Vector V>
    [[nodiscard]] __host__ __device__ inline auto hadamard(V&& x, M&& m) noexcept requires(CUDA<M> && CUDA<V>) {
        return hadamard(m, x);
    }

    template<Matrix T1, Matrix T2>
    [[nodiscard]] inline auto hadamard(T1&& m1, T2&& m2) noexcept requires(CUDA<T1> && CUDA<T2>) {
        return device_obj<MatrixExpr<ExprType::Mul, T1&&, T2&&>>(std::forward<T1>(m1), std::forward<T2>(m2));
    }
}
