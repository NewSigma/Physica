/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.h"

namespace Physica {
    template<Matrix M>
    class device_obj<MatrixExpr<ExprID::Cos, M>>
            : public device_obj<UnaryMatrixExpr<ExprID::Cos, M>> {
        using Base = device_obj<UnaryMatrixExpr<ExprID::Cos, M>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const;

        using Base::reverse;
        void reverse(const Matrix auto& grad) const noexcept;

        [[nodiscard]] auto values(this auto&& self) noexcept;
    };

    template<Matrix M>
    __device__ auto device_obj<MatrixExpr<ExprID::Cos, M>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (isReverseDiff())
            return Base::calc_value(row, col, block);
        else
            return cos(Base::getExpr().calc(row, col, block));
    }

    template<Matrix M>
    void device_obj<MatrixExpr<ExprID::Cos, M>>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& expr = Base::getExpr();
        expr.reverse(-hadamard(grad, sin_elem(expr.values())));
    }

    template<Matrix M>
    auto device_obj<MatrixExpr<ExprID::Cos, M>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return cos_elem(std::forward<Self>(self).getExpr().values());
    }

    template<Matrix M>
    [[nodiscard, gnu::always_inline]] auto cos_elem(M&& m) noexcept requires(DeviceObj<M>) {
        return device_obj<MatrixExpr<ExprID::Cos, remove_device_obj_t<M&&>>>(std::forward<M>(m));
    }
}
