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

#include "../VectorExpr.cuh"

namespace Physica {
    template<Vector T>
    class device_obj<VectorExpr<ExprType::Softmax, T>>
            : public device_obj<UnitaryVectorExpr<ExprType::Softmax, T>> {
        using Base = device_obj<UnitaryVectorExpr<ExprType::Softmax, T>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        template<Vector V, class Executor = SeqExecutor>
        __host__ __device__ inline void assign(V& v) const;

        [[nodiscard]] __device__ ScalarType calc(size_t i) const { return Base::getExpr().softmax(i); }
        [[nodiscard]] __device__ ValueType calc_value(size_t i) const { return Base::getExpr().values().softmax(i); }
        [[nodiscard]] __device__ ScalarType calc(size_t i, ScalarType lnsumexp) const;
        [[nodiscard]] __device__ ValueType calc_value(size_t i, ValueType lnsumexp) const;
    };

    template<Vector T>
    __device__ auto device_obj<VectorExpr<ExprType::Softmax, T>>::calc(size_t i, ScalarType lnsumexp) const -> ScalarType {
        return exp(Base::getExpr().calc(i) - lnsumexp);
    }

    template<Vector T>
    __device__ auto device_obj<VectorExpr<ExprType::Softmax, T>>::calc_value(size_t i, ValueType lnsumexp) const -> ValueType {
        return exp(Base::getExpr().calc_value(i) - lnsumexp);
    }

    template<Vector T>
    template<Vector V, class Executor>
    __host__ __device__ inline void device_obj<VectorExpr<ExprType::Softmax, T>>::assign(V& v) const {
        if constexpr (IsHost())
            Base::assign(v);
        else {
            const auto& expr = Base::getExpr();
            const ScalarType factor = expr.lnSumExp();
            v = exp(expr - factor);
        }
    }

    template<Vector T>
    [[nodiscard]] __host__ __device__ inline auto softmax(T&& v) noexcept requires(CUDA<T>) {
        return device_obj<VectorExpr<ExprType::Softmax, T&&>>(std::forward<T>(v));
    }
}
