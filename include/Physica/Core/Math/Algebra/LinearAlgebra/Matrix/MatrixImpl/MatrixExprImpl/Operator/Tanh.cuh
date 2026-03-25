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
    class device_obj<MatrixExpr<ExprID::Tanh, M>>
            : public device_obj<UnitaryMatrixExpr<ExprID::Tanh, M>> {
        using Base = device_obj<UnitaryMatrixExpr<ExprID::Tanh, M>>;
    public:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const {
            if constexpr (isReverseDiff)
                return calc_value(row, col);
            else
                return tanh(Base::getExpr().calc(row, col));
        }

        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const {
            return tanh(Base::getExpr().calc_value(row, col));
        }

        void reverse(const Matrix auto& y, const Matrix auto& grad) const noexcept;
    };

    template<Matrix M>
    void device_obj<MatrixExpr<ExprID::Tanh, M>>::reverse(const Matrix auto& y, const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& expr = Base::getExpr();
        expr.reverse(hadamard((Tv(1) - square_elem(y)), grad));
    }

    template<Matrix M>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto tanh_elem(M&& m) noexcept requires(DeviceObj<M>) {
        return device_obj<MatrixExpr<ExprID::Tanh, remove_device_obj_t<M&&>>>(std::forward<M>(m));
    }
}
