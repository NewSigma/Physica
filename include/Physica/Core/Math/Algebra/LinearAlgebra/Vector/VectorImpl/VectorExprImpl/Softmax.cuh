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
    template<Vector V>
    class device_obj<VectorExpr<ExprType::Softmax, V>>
            : public device_obj<UnitaryVectorExpr<ExprType::Softmax, V>> {
        using Base = device_obj<UnitaryVectorExpr<ExprType::Softmax, V>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        __host__ __device__ inline void assign(Vector auto& v) const;

        [[nodiscard]] __device__ T calc(size_t i) const { return Base::getExpr().softmax(i); }
        [[nodiscard]] __device__ Tv calc_value(size_t i) const { return Base::getExpr().values().softmax(i); }
        [[nodiscard]] __device__ T calc(size_t i, T lnsumexp) const;
        [[nodiscard]] __device__ Tv calc_value(size_t i, Tv lnsumexp) const;
    };

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprType::Softmax, V>>::calc(size_t i, T lnsumexp) const -> T {
        return exp(Base::getExpr().calc(i) - lnsumexp);
    }

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprType::Softmax, V>>::calc_value(size_t i, Tv lnsumexp) const -> Tv {
        return exp(Base::getExpr().calc_value(i) - lnsumexp);
    }

    template<Vector V>
    template<ExecutePolicy P>
    __host__ __device__ inline void device_obj<VectorExpr<ExprType::Softmax, V>>::assign(Vector auto& v) const {
        if constexpr (IsHost())
            Base::assign(v);
        else {
            const auto& expr = Base::getExpr();
            const T factor = expr.lnSumExp();
            v = exp(expr - factor);
        }
    }

    template<Vector V>
    [[nodiscard]] __host__ __device__ inline auto softmax(V&& v) noexcept requires(CUDA<V>) {
        return device_obj<VectorExpr<ExprType::Softmax, V&&>>(std::forward<V>(v));
    }
}
