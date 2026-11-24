/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.cuh"

namespace Physica {
    template<Vector V>
    class device_obj<VectorExpr<ExprID::Unit, V>>
            : public device_obj<UnaryVectorExpr<ExprID::Unit, V>> {
        using Base = device_obj<UnaryVectorExpr<ExprID::Unit, V>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t i, instanceof_x<ThreadBlock> auto block) const;
    };

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprID::Unit, V>>::calc(size_t i, instanceof_x<ThreadBlock> auto block) const -> T {
        return unit(Base::getExpr().calc(i, block));
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto unit(V&& v) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Unit, remove_device_obj_t<V&&>>>(std::forward<V>(v));
    }
}
