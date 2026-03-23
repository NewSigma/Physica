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
    template<class V, class U>
    class device_obj<VectorExpr<ExprID::Div, V, U>>
            : public device_obj<BinaryVectorExpr<ExprID::Div, V, U>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Div, V, U>>;
    protected:
        using typename Base::T;
    public:
        device_obj(V lhs, U rhs) : Base(std::forward<V>(lhs), std::forward<U>(rhs)) {
            if constexpr (Vector<V>)
                assert(!Base::getRHS().isSubNormal() && "[Error]: Division overflow");
        }
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t index) const {
            if constexpr (Vector<V>)
                return Base::getLHS().calc(index) / Base::getRHS();
            else
                return Base::getLHS() / Base::getRHS().calc(index);
        }
    };

    template<Vector V1, Vector V2>
    class device_obj<VectorExpr<ExprID::Div, V1, V2>>
            : public device_obj<BinaryVectorExpr<ExprID::Div, V1, V2>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Div, V1, V2>>;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t index) const {
            assert(!Base::getRHS().calc(index).isSubNormal() && "[Error]: Division overflow");
            return Base::getLHS().calc(index) / Base::getRHS().calc(index);
        }
    };

    template<Vector V, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator/(V&& v, U&& x) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Div, remove_device_obj_t<V&&>, U&&>>(v, x);
    }

    template<Vector V, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator/(U&& x, V&& v) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Div, U&&, remove_device_obj_t<V&&>>>(x, v);
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto divide(V1&& v1, V2&& v2) noexcept requires(DeviceObj<V1> && DeviceObj<V2>) {
        return device_obj<VectorExpr<ExprID::Div, remove_device_obj_t<V1&&>, remove_device_obj_t<V2&&>>>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
