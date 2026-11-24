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
    class device_obj<VectorExpr<ExprID::Minus, V>>
            : public device_obj<UnaryVectorExpr<ExprID::Minus, V>> {
        using Base = device_obj<UnaryVectorExpr<ExprID::Minus, V>>;
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
    __device__ auto device_obj<VectorExpr<ExprID::Minus, V>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (isReverseDiff())
            return Base::calc_value(index, block);
        else
            return -Base::getExpr().calc(index, block);
    }

    template<Vector V>
    auto device_obj<VectorExpr<ExprID::Minus, V>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return -std::forward<Self>(self).getExpr().values();
    }
}
