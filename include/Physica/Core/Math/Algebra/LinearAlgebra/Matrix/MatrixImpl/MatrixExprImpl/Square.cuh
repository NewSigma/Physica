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

#include "../MatrixExpr.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<MatrixExpr<ExprID::Square, M>>
            : public device_obj<UnitaryMatrixExpr<ExprID::Square, M>> {
        using Base = device_obj<UnitaryMatrixExpr<ExprID::Square, M>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;
    };

    template<Matrix M>
    __device__ auto device_obj<MatrixExpr<ExprID::Square, M>>::calc(size_t row, size_t col) const -> T {
        if constexpr (isReverseDiff)
            return calc_value(row, col);
        else
            return square(Base::getExpr().calc(row, col));
    }

    template<Matrix M>
    __device__ auto device_obj<MatrixExpr<ExprID::Square, M>>::calc_value(size_t row, size_t col) const -> Tv {
        return square(Base::getExpr().calc_value(row, col));
    }

    template<Matrix M>
    void device_obj<MatrixExpr<ExprID::Square, M>>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& expr = Base::getExpr();
        expr.reverse(Tv(2) * hadamard(expr.values(), grad));
    }

    template<Matrix M>
    [[nodiscard]] __host__ __device__ auto square_elem(M&& m) noexcept requires(CUDA<M>) {
        return device_obj<MatrixExpr<ExprID::Square, remove_device_obj_t<M&&>>>(std::forward<M>(m));
    }
}
