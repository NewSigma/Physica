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
    class device_obj<MatrixExpr<ExprID::Div, T, U>>
            : public device_obj<BinaryMatrixExpr<ExprID::Div, T, U>> {
        using Base = BinaryMatrixExpr<ExprID::Div, T, U>;
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
    class device_obj<MatrixExpr<ExprID::Div, T, U>>
            : public device_obj<BinaryMatrixExpr<ExprID::Div, T, U>> {
        using Base = device_obj<BinaryMatrixExpr<ExprID::Div, T, U>>;
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

    template<Matrix M, Matrix U>
    class device_obj<MatrixExpr<ExprID::Div, M, U>>
            : public device_obj<BinaryMatrixExpr<ExprID::Div, M, U>> {
        using Base = device_obj<BinaryMatrixExpr<ExprID::Div, M, U>>;
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

    template<Matrix M, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator/(M&& m, U&& x) noexcept requires(CUDA<M>) {
        return device_obj<MatrixExpr<ExprID::Div, M&&, U&&>>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator/(U&& x, M&& m) noexcept requires(CUDA<M>) {
        return device_obj<MatrixExpr<ExprID::Div, U&&, M&&>>(std::forward<U>(x), std::forward<M>(m));
    }

    template<Matrix M, Vector U>
    [[nodiscard]] __host__ __device__ auto divide(M&& m, U&& x) noexcept requires(CUDA<M> && CUDA<U>) {
        return device_obj<MatrixExpr<ExprID::Div, M&&, U&&>>(std::forward<M>(m), std::forward<U>(x));
    }

    template<Matrix M, Vector U>
    [[nodiscard]] __host__ __device__ auto divide(U&& x, M&& m) noexcept requires(CUDA<M> && CUDA<U>) {
        return device_obj<MatrixExpr<ExprID::Div, U&&, M&&>>(std::forward<U>(x), std::forward<M>(m));
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard]] __host__ __device__ auto divide(M1&& m1, M2&& m2) noexcept requires(CUDA<M1> && CUDA<M2>) {
        return device_obj<MatrixExpr<ExprID::Div, M1&&, M2&&>>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}
