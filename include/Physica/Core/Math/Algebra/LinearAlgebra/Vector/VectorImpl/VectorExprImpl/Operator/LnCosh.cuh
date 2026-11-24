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
    class device_obj<VectorExpr<ExprID::LnCosh, V>>
            : public device_obj<UnaryVectorExpr<ExprID::LnCosh, V>> {
        using Base = device_obj<UnaryVectorExpr<ExprID::LnCosh, V>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;

        [[nodiscard]] auto values(this auto&& self) noexcept;
    };

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprID::LnCosh, V>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (isReverseDiff())
            return Base::calc_value(index, block);
        else
            return lncosh(Base::getExpr().calc(index, block));
    }

    template<Vector V>
    auto device_obj<VectorExpr<ExprID::LnCosh, V>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return lncosh(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto lncosh(V&& v) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::LnCosh, remove_device_obj_t<V&&>>>(std::forward<V>(v));
    }
}
