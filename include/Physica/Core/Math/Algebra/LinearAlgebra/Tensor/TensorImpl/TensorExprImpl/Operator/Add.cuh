/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/TensorImpl/TensorExpr.cuh"

namespace Physica {
    template<Tensor X1, Tensor X2>
    class device_obj<TensorExpr<ExprID::Add, X1, X2>>
            : public device_obj<BinaryTensorExpr<ExprID::Add, X1, X2>> {
        using Base = device_obj<BinaryTensorExpr<ExprID::Add, X1, X2>>;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] __device__ T calc(const typename Base::IndexType& indices) const {
            return Base::getLHS().calc(indices) + Base::getRHS().calc(indices);
        }
    };

    template<Tensor X, Tensor Y>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator+(X&& x, Y&& y) noexcept requires(DeviceObj<X> && DeviceObj<Y>) {
        if constexpr (!canonicalized(x, y))
            return std::forward<Y>(y) + std::forward<X>(x);
        else
            return device_obj<TensorExpr<ExprID::Add, remove_device_obj_t<X&&>, remove_device_obj_t<Y&&>>>(std::forward<X>(x), std::forward<Y>(y));
    }
}
