/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.cuh"

namespace Physica {
    template<Vector V>
    class device_obj<VectorExpr<ExprID::Cos, V>>
            : public device_obj<UnaryVectorExpr<ExprID::Cos, V>> {
        using Base = device_obj<UnaryVectorExpr<ExprID::Cos, V>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Getters */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;

        [[nodiscard]] auto values(this auto&& self) noexcept;
    };

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprID::Cos, V>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (isReverseDiff())
            return Base::calc_value(index, block);
        else
            return cos(Base::getExpr().calc(index, block));
    }

    template<Vector V>
    auto device_obj<VectorExpr<ExprID::Cos, V>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return cos(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto cos(V&& v) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Cos, remove_device_obj_t<V&&>>>(std::forward<V>(v));
    }
}
