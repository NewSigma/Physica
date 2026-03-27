/*
 * Copyright 2025-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.cuh"

namespace Physica {
    template<Matrix M1, Matrix M2>
    class device_obj<MatrixExpr<ExprID::Sub, M1, M2>>
            : public device_obj<BinaryMatrixExpr<ExprID::Sub, M1, M2>> {
        using host_obj = MatrixExpr<ExprID::Sub, M1, M2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<BinaryMatrixExpr<ExprID::Sub, M1, M2>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const;
    };

    template<Matrix M1, Matrix M2>
    __device__ auto device_obj<MatrixExpr<ExprID::Sub, M1, M2>>::calc(size_t row, size_t col) const -> T {
        return Base::getLHS().calc(row, col) - Base::getRHS().calc(row, col);
    }

    template<Matrix M1, Matrix M2>
    __device__ auto device_obj<MatrixExpr<ExprID::Sub, M1, M2>>::calc_value(size_t row, size_t col) const -> Tv {
        return Base::getLHS().calc_value(row, col) - Base::getRHS().calc_value(row, col);
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator-(M&& m, U&& x) noexcept requires(DeviceObj<M>) {
        return std::forward<M>(m) + (-std::forward<U>(x));
    }

    template<Matrix M, Scalar U>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator-(U&& x, M&& m) noexcept requires(DeviceObj<M>) {
        return std::forward<U>(x) + (-std::forward<M>(m));
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator-(M1&& m1, M2&& m2) noexcept requires(DeviceObj<M1> && DeviceObj<M2>){
        return device_obj<MatrixExpr<ExprID::Sub, remove_device_obj_t<M1&&>, remove_device_obj_t<M2&&>>>(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}
