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

#include "../TensorExpr.h"

namespace Physica {
    template<Tensor T1, Tensor T2>
    class TensorExpr<ExprType::Add, T1, T2>
            : public BinaryTensorExpr<ExprType::Add, T1, T2> {
        using Base = BinaryTensorExpr<ExprType::Add, T1, T2>;
    public:
        using typename Base::ScalarType;
        using typename Base::IndexArray;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] ScalarType calc(const IndexArray& indices) const {
            return Base::getLHS().calc(indices) + Base::getRHS().calc(indices);
        }
    };

    template<Tensor T1, Tensor T2>
    [[nodiscard]] inline auto operator+(const T1& x, const T2& y) noexcept {
        return TensorExpr<ExprType::Add, T1, T2>(x, y);
    }
}
