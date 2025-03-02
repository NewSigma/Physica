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
    template<class T, class U>
    class device_obj<MatrixExpr<ExprType::Sub, T, U>> : public device_obj<BinaryMatrixExpr<ExprType::Sub, T, U>> {
        static_assert(Scalar<T> || Scalar<U>, "[Error]: Either types should be Scalar");

        using Base = device_obj<BinaryMatrixExpr<ExprType::Sub, T, U>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc(row, col) - Base::getRHS();
            else
                return Base::getLHS() - Base::getRHS().calc(row, col);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            if constexpr (Matrix<T>)
                return Base::getLHS().calc_value(row, col) - Base::getRHS().value();
            else
                return Base::getLHS().value() - Base::getRHS().calc_value(row, col);
        }
    };

    template<Matrix T1, Matrix T2>
    class device_obj<MatrixExpr<ExprType::Sub, T1, T2>>
            : public device_obj<BinaryMatrixExpr<ExprType::Sub, T1, T2>> {
        using host_obj = MatrixExpr<ExprType::Sub, T1, T2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<BinaryMatrixExpr<ExprType::Sub, T1, T2>>;
        constexpr static bool IsSymm = MatrixOption::isSymmMatrix<T1>() && MatrixOption::isSymmMatrix<T2>();
        constexpr static bool IsHermite = MatrixOption::isHermiteMatrix<T1>() && MatrixOption::isHermiteMatrix<T2>();
        using TransposeRtnTy = std::conditional<IsSymm, const This&, device_obj<Transpose<host_obj>>>::type;
        using HermiteRtnTy = std::conditional<IsHermite, const This&, device_obj<Hermite<host_obj>>>::type;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const {
            return Base::getLHS().calc(row, col) - Base::getRHS().calc(row, col);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            return Base::getLHS().calc_value(row, col) - Base::getRHS().calc_value(row, col);
        }

        [[nodiscard]] TransposeRtnTy transpose() const noexcept { return TransposeRtnTy(*this); }
        [[nodiscard]] HermiteRtnTy hermite() const noexcept { return HermiteRtnTy(*this); }
    };

    template<Matrix T, Scalar U>
    [[nodiscard]] __host__ __device__ inline auto operator-(T&& m, U&& x) noexcept requires(CUDA<T>) {
        return device_obj<MatrixExpr<ExprType::Sub, T&&, U&&>>(std::forward<T>(m), std::forward<U>(x));
    }

    template<Matrix T, Scalar U>
    [[nodiscard]] __host__ __device__ inline auto operator-(U&& x, T&& m) noexcept requires(CUDA<T>) {
        return device_obj<MatrixExpr<ExprType::Sub, U&&, T&&>>(std::forward<U>(x), std::forward<T>(m));
    }

    template<Matrix T1, Matrix T2>
    [[nodiscard]] __host__ __device__ inline auto operator-(T1&& m1, T2&& m2) noexcept requires(CUDA<T1> && CUDA<T2>){
        return device_obj<MatrixExpr<ExprType::Sub, T1&&, T2&&>>(std::forward<T1>(m1), std::forward<T2>(m2));
    }
}
