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

#include "../VectorExpr.cuh"

namespace Physica {
    template<Vector T, Scalar U>
    class device_obj<VectorExpr<ExprID::Sub, T, U>>
            : public device_obj<BinaryVectorExpr<ExprID::Sub, T, U>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Sub, T, U>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return Base::getLHS().calc(index) - Base::getRHS();
        }

    };

    template<Vector T1, Vector T2>
    class device_obj<VectorExpr<ExprID::Sub, T1, T2>>
            : public device_obj<BinaryVectorExpr<ExprID::Sub, T1, T2>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Sub, T1, T2>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return Base::getLHS().calc(index) - Base::getRHS().calc(index);
        }
    };

    template<Vector T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator-(T&& v, U&& x) noexcept requires(CUDA<T>) {
        return device_obj<VectorExpr<ExprID::Sub, T&&, U&&>>(std::forward<T>(v), std::forward<U>(x));
    }

    template<Vector T1, Vector T2>
    [[nodiscard]] __host__ __device__ auto operator-(T1&& v1, T2&& v2) noexcept requires(CUDA<T1> && CUDA<T2>) {
        return device_obj<VectorExpr<ExprID::Sub, T1&&, T2&&>>(std::forward<T1>(v1), std::forward<T2>(v2));
    }
}
