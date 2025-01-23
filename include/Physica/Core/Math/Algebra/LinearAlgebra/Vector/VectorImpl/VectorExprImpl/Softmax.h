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

#include "../VectorExpr.h"

namespace Physica::Core {
    template<Vector T>
    class VectorExpr<ExprType::Softmax, T>
            : public UnitaryVectorExpr<ExprType::Softmax, T> {
        using Base = UnitaryVectorExpr<ExprType::Softmax, T>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t i) const { return Base::getExpr().softmax(i); }

        [[nodiscard]] ValueType calc_value(size_t i) const { return Base::getExpr().values().softmax(i); }

        template<Vector U, Vector V>
        void reverse(const U& y, const V& grad_) const noexcept requires(isReverseDiff) ;
    };

    template<Vector T>
    template<Vector U, Vector V>
    void VectorExpr<ExprType::Softmax, T>::reverse(const U& y, const V& grad_) const noexcept requires(isReverseDiff) {
        const auto& grad = grad_.values();
        Base::getExpr().reverse(hadamard(y, grad) - (y * grad) * y);
    }

    template<Vector T>
    [[nodiscard]] inline auto softmax(const T& v) noexcept {
        return VectorExpr<ExprType::Softmax, T>(v);
    }
}
