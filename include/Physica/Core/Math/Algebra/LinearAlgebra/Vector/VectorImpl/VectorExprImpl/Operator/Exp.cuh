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
    class device_obj<VectorExpr<ExprID::Exp, V>>
            : public device_obj<UnaryVectorExpr<ExprID::Exp, V>> {
        using Base = device_obj<UnaryVectorExpr<ExprID::Exp, V>>;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;
    };

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprID::Exp, V>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        return exp(Base::getExpr().calc(index, block));
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto exp(V&& v) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Exp, remove_device_obj_t<V&&>>>(std::forward<V>(v));
    }
}
