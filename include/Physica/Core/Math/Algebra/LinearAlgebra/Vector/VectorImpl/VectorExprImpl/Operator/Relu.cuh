/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.cuh"

namespace Physica {
    template<Vector V>
    class device_obj<VectorExpr<ExprID::Relu, V>>
            : public device_obj<UnaryVectorExpr<ExprID::Relu, V>> {
        using Base = device_obj<UnaryVectorExpr<ExprID::Relu, V>>;
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

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&& self) noexcept;
    };

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprID::Relu, V>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (isReverseDiff())
            return Base::calc_value(index, block);
        else
            return relu(Base::getExpr().calc(index, block));
    }

    template<Vector V>
    void device_obj<VectorExpr<ExprID::Relu, V>>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& expr = Base::getExpr();
        expr.reverse(hadamard(relu(unit(expr.values())), grad.values()));
    }

    template<Vector V>
    auto device_obj<VectorExpr<ExprID::Relu, V>>::values(this auto&& self) noexcept {
        return relu(std::forward<decltype(self)>(self).getExpr().values());
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto relu(V&& v) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Relu, remove_device_obj_t<V&&>>>(std::forward<V>(v));
    }
}
