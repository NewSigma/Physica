/*
 * Copyright 2025 Weibo He.
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

#include "../MatrixExpr.h"

namespace Physica {
    template<Matrix M>
    class device_obj<MatrixExpr<ExprType::Minus, M>>
            : public device_obj<UnitaryMatrixExpr<ExprType::Minus, M>> {
        using Base = device_obj<UnitaryMatrixExpr<ExprType::Minus, M>>;
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

        using Base::reverse;
        void reverse(const Matrix auto& grad) const noexcept;
    };

    template<Matrix M>
    __device__ auto device_obj<MatrixExpr<ExprType::Minus, M>>::calc(size_t row, size_t col) const -> T {
        if constexpr (isReverseDiff)
            return calc_value(row, col);
        else
            return -Base::getExpr().calc(row, col);
    }

    template<Matrix M>
    __device__ auto device_obj<MatrixExpr<ExprType::Minus, M>>::calc_value(size_t row, size_t col) const -> Tv {
        return -Base::getExpr().calc_value(row, col);
    }

    template<Matrix M>
    void device_obj<MatrixExpr<ExprType::Minus, M>>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff);
        Base::getExpr().reverse(-grad);
    }

    template<Matrix M>
    [[nodiscard]] auto operator-(M&& m) noexcept requires(CUDA<M>) {
        return device_obj<MatrixExpr<ExprType::Minus, M&&>>(std::forward<M>(m));
    }
}
