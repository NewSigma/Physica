/*
 * Copyright 2024-2025 Weibo He.
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

namespace Physica {
    template<Vector T>
    class VectorExpr<ExprType::Relu, T> : public UnitaryVectorExpr<ExprType::Relu, T> {
        using This = VectorExpr<ExprType::Relu, T>;
        using Base = UnitaryVectorExpr<ExprType::Relu, T>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t index) const {
            return relu(Base::getExpr().calc(index));
        }

        [[nodiscard]] ValueType calc_value(size_t index) const {
            return relu(Base::getExpr().calc_value(index));
        }

        template<Vector V>
        void reverse(const V& grad_) const noexcept requires(isReverseDiff) {
            const auto& grad = grad_.values();
            for (size_t i = 0; i < Base::getLength(); ++i) {
                auto v = Base::getExpr().calc(i);
                v.reverse(v.isPositive() ? grad.calc(i) : ValueType(0));
            }
        }
    };

    template<Vector T>
    [[nodiscard]] inline auto relu(T&& v) noexcept requires(!CUDA<T>) {
        return VectorExpr<ExprType::Relu, T&&>(std::forward<T>(v));
    }
}
