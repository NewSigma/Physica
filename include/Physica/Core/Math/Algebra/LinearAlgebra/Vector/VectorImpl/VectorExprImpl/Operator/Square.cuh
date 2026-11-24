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
    class device_obj<VectorExpr<ExprID::Square, V>>
            : public device_obj<UnaryVectorExpr<ExprID::Square, V>> {
        using Base = device_obj<UnaryVectorExpr<ExprID::Square, V>>;
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

        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&& self) noexcept;
    };

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprID::Square, V>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (Base::isReverseDiff())
            return Base::calc_value(index, block);
        else
            return square(Base::getExpr().calc(index, block));
    }

    template<Vector V>
    void device_obj<VectorExpr<ExprID::Square, V>>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& expr = Base::getExpr();
        expr.reverse(expr.values() * (Tv(2) * grad));
    }

    template<Vector V>
    auto device_obj<VectorExpr<ExprID::Square, V>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return square(std::forward<Self>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto square(V&& v) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Square, remove_device_obj_t<V&&>>>(std::forward<V>(v));
    }
}
