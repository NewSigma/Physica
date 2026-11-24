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
    class device_obj<MatrixExpr<ExprID::Unit, M>>
            : public device_obj<UnaryMatrixExpr<ExprID::Unit, M>> {
        using Base = device_obj<UnaryMatrixExpr<ExprID::Unit, M>>;
    public:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const;

        [[nodiscard]] auto values(this auto&& self) noexcept;
    };

    template<Matrix M>
    __device__ auto device_obj<MatrixExpr<ExprID::Unit, M>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (isReverseDiff())
            return Base::calc_value(row, col, block);
        else
            return unit(Base::getExpr().calc(row, col, block));
    }

    template<Matrix M>
    auto device_obj<MatrixExpr<ExprID::Unit, M>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return unit_elem(std::forward<Self>(self).getExpr().values());
    }

    template<Matrix M>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto unit_elem(M&& m) noexcept requires(DeviceObj<M>) {
        return device_obj<MatrixExpr<ExprID::Unit, remove_device_obj_t<M&&>>>(std::forward<M>(m));
    }
}
