/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<MatrixExpr<ExprID::Sinh, M>>
            : public device_obj<UnaryMatrixExpr<ExprID::Sinh, M>> {
        using Base = device_obj<UnaryMatrixExpr<ExprID::Sinh, M>>;
    public:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const;

        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] auto values(this auto&& self) noexcept;
    };

    template<Matrix M>
    __device__ auto device_obj<MatrixExpr<ExprID::Sinh, M>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (isReverseDiff())
            return Base::calc_value(row, col, block);
        else
            return sinh(Base::getExpr().calc(row, col, block));
    }

    template<Matrix M>
    void device_obj<MatrixExpr<ExprID::Sinh, M>>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& expr = Base::getExpr();
        expr.reverse(hadamard(cosh_elem(expr.values()), grad));
    }

    template<Matrix M>
    auto device_obj<MatrixExpr<ExprID::Sinh, M>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return sinh_elem(std::forward<Self>(self).getExpr().values());
    }

    template<Matrix M>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto sinh_elem(M&& m) noexcept requires(DeviceObj<M>) {
        return device_obj<MatrixExpr<ExprID::Sinh, remove_device_obj_t<M&&>>>(std::forward<M>(m));
    }
}
