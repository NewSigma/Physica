/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<MatrixExpr<ExprID::Sqrt, M>>
            : public device_obj<UnaryMatrixExpr<ExprID::Sqrt, M>> {
        using Base = device_obj<UnaryMatrixExpr<ExprID::Sqrt, M>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const;
        [[nodiscard]] auto values(this auto&& self) noexcept;
    };

    template<Matrix M>
    __device__ auto device_obj<MatrixExpr<ExprID::Sqrt, M>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (Base::isReverseDiff())
            return Base::calc_value(row, col, block);
        else
            return sqrt(Base::getExpr().calc(row, col, block));
    }

    template<Matrix M>
    auto device_obj<MatrixExpr<ExprID::Sqrt, M>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return sqrt_elem(std::forward<Self>(self).getExpr().values());
    }

    template<Matrix M>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto sqrt_elem(M&& m) noexcept requires(DeviceObj<M>) {
        return device_obj<MatrixExpr<ExprID::Sqrt, remove_device_obj_t<M&&>>>(std::forward<M>(m));
    }
}
