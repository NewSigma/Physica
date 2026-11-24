/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
