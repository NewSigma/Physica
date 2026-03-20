/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.cuh"

namespace Physica {
    template<Vector V>
    class device_obj<VectorExpr<ExprID::Reciprocal, V>>
            : public device_obj<UnitaryVectorExpr<ExprID::Reciprocal, V>> {
        using Base = device_obj<UnitaryVectorExpr<ExprID::Reciprocal, V>>;
    public:
        using typename Base::ScalarType;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            return reciprocal(Base::getExpr().calc(index));
        }

        template<Packet Pack>
        [[nodiscard]] __device__ Pack packet(size_t index) const noexcept {
            return Pack(1) / Base::getExpr().template packet<Pack>(index);
        }

        template<Packet Pack>
        [[nodiscard]] __device__ Pack packet(size_t index, size_t count) const noexcept {
            return Pack(1) / Base::getExpr().template packet<Pack>(index, count);
        }
    };

    template<Vector V>
    [[nodiscard]] __host__ __device__ auto reciprocal(V&& v) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Reciprocal, remove_device_obj_t<V&&>>>(std::forward<V>(v));
    }
}
