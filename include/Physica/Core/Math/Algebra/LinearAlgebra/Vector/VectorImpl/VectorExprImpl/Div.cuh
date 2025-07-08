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
    template<class T, class U>
    class device_obj<VectorExpr<ExprType::Div, T, U>>
            : public device_obj<BinaryVectorExpr<ExprType::Div, T, U>> {
        using Base = device_obj<BinaryVectorExpr<ExprType::Div, T, U>>;
    public:
        using typename Base::ScalarType;
    public:
        device_obj(T lhs, U rhs) : Base(std::forward<T>(lhs), std::forward<U>(rhs)) {
            if constexpr (Vector<T>)
                assert(!Base::getRHS().isZero() && "[Error]: Divide by zero");
        }
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            if constexpr (Vector<T>)
                return Base::getLHS().calc(index) / Base::getRHS();
            else
                return Base::getLHS() / Base::getRHS().calc(index);
        }
    };

    template<Vector T1, Vector T2>
    class device_obj<VectorExpr<ExprType::Div, T1, T2>>
            : public device_obj<BinaryVectorExpr<ExprType::Div, T1, T2>> {
        using Base = device_obj<BinaryVectorExpr<ExprType::Div, T1, T2>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            assert(!Base::getRHS().calc(index).isZero() && "[Error]: Divide by zero");
            return Base::getLHS().calc(index) / Base::getRHS().calc(index);
        }
    };

    template<Vector T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator/(T&& v, U&& x) noexcept requires(CUDA<T>) {
        return device_obj<VectorExpr<ExprType::Div, T&&, U&&>>(v, x);
    }

    template<Vector T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator/(U&& x, T&& v) noexcept requires(CUDA<T>) {
        return device_obj<VectorExpr<ExprType::Div, U&&, T&&>>(x, v);
    }

    template<Vector T1, Vector T2>
    [[nodiscard]] auto divide(T1&& v1, T2&& v2) noexcept requires(CUDA<T1> && CUDA<T2>) {
        return device_obj<VectorExpr<ExprType::Div, T1&&, T2&&>>(std::forward<T1>(v1), std::forward<T2>(v2));
    }
}
