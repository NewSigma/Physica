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
    template<Matrix T, Scalar U>
    class device_obj<MatrixExpr<ExprType::Add, T, U>>
            : public device_obj<BinaryMatrixExpr<ExprType::Add, T, U>> {
        using Base = device_obj<BinaryMatrixExpr<ExprType::Add, T, U>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const {
            return Base::getLHS().calc(row, col) + Base::getRHS();
        }

        __host__ __device__ auto values() const noexcept { return Base::getLHS().values() + Base::getRHS().values(); }
    };

    template<Matrix T, Vector U>
    class device_obj<MatrixExpr<ExprType::Add, T, U>>
            : public device_obj<BinaryMatrixExpr<ExprType::Add, T, U>> {
        using Base = device_obj<BinaryMatrixExpr<ExprType::Add, T, U>>;
    public:
        using typename Base::ScalarType;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        template<Matrix M>
        __host__ __device__ void assign(M& target) const {
            constexpr bool FastAssign = Traits<std::remove_cvref_t<T>>::FastAssign;
            if constexpr (FastAssign) {
                Base::getLHS().template assign<M>(target);
                target += Base::getRHS();
            }
            else
                Base::template assign<M>(target);
        }

        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const {
            if constexpr (isReverseDiff)
                return calc_value(row, col);
            else
                return Base::getLHS().calc(row, col) + Base::getRHS().calc(row);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) + Base::getRHS().calc_value(row);
        }

        template<Matrix M1>
        void reverse(const M1& grad) const noexcept requires(isReverseDiff) {
            if constexpr (ReverseDiff<T>)
                Base::getLHS().reverse(grad);
            if constexpr (ReverseDiff<U>)
                Base::getRHS().reverse(grad.sum_cols());
        }

        __host__ __device__ auto values() const noexcept { return Base::getLHS().values() + Base::getRHS().values(); }
    };

    template<Matrix T1, Matrix T2>
    class device_obj<MatrixExpr<ExprType::Add, T1, T2>>
            : public device_obj<BinaryMatrixExpr<ExprType::Add, T1, T2>> {
        using Base = device_obj<BinaryMatrixExpr<ExprType::Add, T1, T2>>;
    public:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const {
            if constexpr (isReverseDiff)
                return calc_value(row, col);
            else
                return Base::getLHS().calc(row, col) + Base::getRHS().calc(row, col);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) + Base::getRHS().calc_value(row, col);
        }

        __host__ __device__ auto values() const noexcept { return Base::getLHS().values() + Base::getRHS().values(); }
    };

    template<Matrix T, Scalar U>
    [[nodiscard]] __host__ __device__ inline auto operator+(T&& m, U&& x) noexcept requires(CUDA<T>) {
        return device_obj<MatrixExpr<ExprType::Add, T&&, U&&>>(std::forward<T>(m), std::forward<U>(x));
    }

    template<Matrix T, Scalar U>
    [[nodiscard]] __host__ __device__ inline auto operator+(U&& x, T&& m) noexcept requires(CUDA<T>) {
        return m + x;
    }

    template<Matrix T, Vector U>
    [[nodiscard]] __host__ __device__ inline auto operator+(T&& m, U&& x) noexcept requires(CUDA<T> && CUDA<U>) {
        return device_obj<MatrixExpr<ExprType::Add, T&&, U&&>>(std::forward<T>(m), std::forward<U>(x));
    }

    template<Matrix T, Vector U>
    [[nodiscard]] __host__ __device__ inline auto operator+(U&& x, T&& m) noexcept requires(CUDA<T> && CUDA<U>) {
        return m + x;
    }

    template<Matrix T1, Matrix T2>
    [[nodiscard]] __host__ __device__ inline auto operator+(T1&& m1, T2&& m2) noexcept requires(CUDA<T1> && CUDA<T2>) {
        return device_obj<MatrixExpr<ExprType::Add, T1&&, T2&&>>(std::forward<T1>(m1), std::forward<T2>(m2));
    }
}
