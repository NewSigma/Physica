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
    template<Vector V, Scalar U>
    class device_obj<VectorExpr<ExprID::Sub, V, U>>
            : public device_obj<BinaryVectorExpr<ExprID::Sub, V, U>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Sub, V, U>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return Base::getLHS().calc(index) - Base::getRHS();
        }

    };

    template<Vector V1, Vector V2>
    class device_obj<VectorExpr<ExprID::Sub, V1, V2>>
            : public device_obj<BinaryVectorExpr<ExprID::Sub, V1, V2>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Sub, V1, V2>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return Base::getLHS().calc(index) - Base::getRHS().calc(index);
        }
    };

    template<Vector V, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator-(V&& v, U&& x) noexcept requires(CUDA<V>) {
        return device_obj<VectorExpr<ExprID::Sub, V&&, U&&>>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] __host__ __device__ auto operator-(V1&& v1, V2&& v2) noexcept requires(CUDA<V1> && CUDA<V2>) {
        return device_obj<VectorExpr<ExprID::Sub, V1&&, V2&&>>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
