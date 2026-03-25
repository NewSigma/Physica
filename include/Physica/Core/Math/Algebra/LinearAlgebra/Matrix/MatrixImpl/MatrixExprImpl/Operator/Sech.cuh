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
    class device_obj<MatrixExpr<ExprID::Sech, M>>
            : public device_obj<UnitaryMatrixExpr<ExprID::Sech, M>> {
        using Base = device_obj<UnitaryMatrixExpr<ExprID::Sech, M>>;
    public:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const;
    };

    template<Matrix M>
    __device__ auto device_obj<MatrixExpr<ExprID::Sech, M>>::calc(size_t row, size_t col) const -> T {
        if constexpr (isReverseDiff)
            return calc_value(row, col);
        else
            return sech(Base::getExpr().calc(row, col));
    }

    template<Matrix M>
    __device__ auto device_obj<MatrixExpr<ExprID::Sech, M>>::calc_value(size_t row, size_t col) const -> Tv {
        return sech(Base::getExpr().calc_value(row, col));
    }

    template<Matrix M>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto sech_elem(M&& m) noexcept requires(DeviceObj<M>) {
        return device_obj<MatrixExpr<ExprID::Sech, remove_device_obj_t<M&&>>>(std::forward<M>(m));
    }
}
