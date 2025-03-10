/*
 * Copyright 2025 Weibo He.
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
    template<class T, class U> requires(Scalar<T> || Scalar<U>)
    class device_obj<MatrixExpr<ExprType::Div, T, U>>
            : public device_obj<BinaryMatrixExpr<ExprType::Div, T, U>> {
        using Base = BinaryMatrixExpr<ExprType::Div, T, U>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::Tv;
    public:
        device_obj(T lhs, U rhs) : Base(std::forward<T>(lhs), std::forward<U>(rhs)) {
            if constexpr (Matrix<T>)
                assert(!Base::getRHS().isZero() && "[Error]: Divide by zero");
        }
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc(row, col) / Base::getRHS();
            else
                return Base::getLHS() / Base::getRHS().calc(row, col);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc_value(row, col) / Base::getRHS().value();
            else
                return Base::getLHS().value() / Base::getRHS().calc_value(row, col);
        }
    };

    template<class T, class U> requires(Vector<T> || Vector<U>)
    class device_obj<MatrixExpr<ExprType::Div, T, U>>
            : public device_obj<BinaryMatrixExpr<ExprType::Div, T, U>> {
        using Base = device_obj<BinaryMatrixExpr<ExprType::Div, T, U>>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc(row, col) / Base::getRHS().calc(row);
            else
                return Base::getLHS().calc(row) / Base::getRHS().calc(row, col);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc_value(row, col) / Base::getRHS().calc_value(row);
            else
                return Base::getLHS().calc_value(row) / Base::getRHS().calc_value(row, col);
        }
    };

    template<Matrix T, Matrix U>
    class device_obj<MatrixExpr<ExprType::Div, T, U>>
            : public device_obj<BinaryMatrixExpr<ExprType::Div, T, U>> {
        using Base = device_obj<BinaryMatrixExpr<ExprType::Div, T, U>>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const {
            assert(!Base::getRHS().calc(row, col).isZero() && "[Error]: Divide by zero");
            return Base::getLHS().calc(row, col) / Base::getRHS().calc(row, col);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            assert(!Base::getRHS().calc_value(row, col).isZero() && "[Error]: Divide by zero");
            return Base::getLHS().calc_value(row, col) / Base::getRHS().calc_value(row, col);
        }
    };

    template<Matrix T, Scalar U>
    [[nodiscard]] __host__ __device__ inline auto operator/(T&& m, U&& x) noexcept requires(CUDA<T>) {
        return device_obj<MatrixExpr<ExprType::Div, T&&, U&&>>(std::forward<T>(m), std::forward<U>(x));
    }

    template<Matrix T, Scalar U>
    [[nodiscard]] __host__ __device__ inline auto operator/(U&& x, T&& m) noexcept requires(CUDA<T>) {
        return device_obj<MatrixExpr<ExprType::Div, U&&, T&&>>(std::forward<U>(x), std::forward<T>(m));
    }

    template<Matrix T, Vector U>
    [[nodiscard]] __host__ __device__ inline auto divide(T&& m, U&& x) noexcept requires(CUDA<T> && CUDA<U>) {
        return device_obj<MatrixExpr<ExprType::Div, T&&, U&&>>(std::forward<T>(m), std::forward<U>(x));
    }

    template<Matrix T, Vector U>
    [[nodiscard]] __host__ __device__ inline auto divide(U&& x, T&& m) noexcept requires(CUDA<T> && CUDA<U>) {
        return device_obj<MatrixExpr<ExprType::Div, U&&, T&&>>(std::forward<U>(x), std::forward<T>(m));
    }

    template<Matrix T1, Matrix T2>
    [[nodiscard]] __host__ __device__ inline auto divide(T1&& m1, T2&& m2) noexcept requires(CUDA<T1> && CUDA<T2>) {
        return device_obj<MatrixExpr<ExprType::Div, T1&&, T2&&>>(std::forward<T1>(m1), std::forward<T2>(m2));
    }
}
