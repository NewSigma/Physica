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
    class device_obj<VectorExpr<ExprID::Exp, V>>
            : public device_obj<UnaryVectorExpr<ExprID::Exp, V>> {
        using Base = device_obj<UnaryVectorExpr<ExprID::Exp, V>>;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;
    };

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprID::Exp, V>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        return exp(Base::getExpr().calc(index, block));
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto exp(V&& v) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Exp, remove_device_obj_t<V&&>>>(std::forward<V>(v));
    }
}
