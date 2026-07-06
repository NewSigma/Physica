/*
 * Copyright 2025-2026 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<MatrixExpr<ExprID::Cosh, M>>
            : public device_obj<UnitaryMatrixExpr<ExprID::Cosh, M>> {
        using Base = device_obj<UnitaryMatrixExpr<ExprID::Cosh, M>>;
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
    __device__ auto device_obj<MatrixExpr<ExprID::Cosh, M>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (isReverseDiff())
            return Base::calc_value(row, col, block);
        else
            return cosh(Base::getExpr().calc(row, col, block));
    }

    template<Matrix M>
    void device_obj<MatrixExpr<ExprID::Cosh, M>>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff());
        const auto& expr = Base::getExpr();
        expr.reverse(hadamard(sinh_elem(expr.values()), grad));
    }

    template<Matrix M>
    auto device_obj<MatrixExpr<ExprID::Cosh, M>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return cosh_elem(std::forward<Self>(self).getExpr().values());
    }

    template<Matrix M>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto cosh_elem(M&& m) noexcept requires(DeviceObj<M>) {
        return device_obj<MatrixExpr<ExprID::Cosh, remove_device_obj_t<M&&>>>(std::forward<M>(m));
    }
}
