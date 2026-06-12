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
    template<class U, class V>
    class device_obj<VectorExpr<ExprID::Div, U, V>>
            : public device_obj<BinaryVectorExpr<ExprID::Div, U, V>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Div, U, V>>;
    protected:
        using typename Base::T;
        using typename Base::Ref1;
        using typename Base::Ref2;
    public:
        __host__ __device__ device_obj(Ref1 lhs, Ref2 rhs);
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;
    };

    template<class U, class V>
    __host__ __device__ device_obj<VectorExpr<ExprID::Div, U, V>>::device_obj(Ref1 lhs, Ref2 rhs) : Base(std::forward<Ref1>(lhs), std::forward<Ref2>(rhs)) {
        if constexpr (Vector<U>)
            assert(!Base::getRHS().isSubNormal() && "[Error]: Division overflow");
    }

    template<class U, class V>
    __device__ auto device_obj<VectorExpr<ExprID::Div, U, V>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        if constexpr (Vector<U>)
            return Base::getLHS().calc(index, block) / Base::getRHS();
        else
            return Base::getLHS() / Base::getRHS().calc(index, block);
    }

    template<Vector V1, Vector V2>
    class device_obj<VectorExpr<ExprID::Div, V1, V2>>
            : public device_obj<BinaryVectorExpr<ExprID::Div, V1, V2>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Div, V1, V2>>;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;
    };

    template<Vector V1, Vector V2>
    __device__ auto device_obj<VectorExpr<ExprID::Div, V1, V2>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        assert(!Base::getRHS().calc(index, block).isSubNormal() && "[Error]: Division overflow");
        return Base::getLHS().calc(index, block) / Base::getRHS().calc(index, block);
    }

    template<Vector V, Scalar U>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator/(V&& v, U&& x) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Div, remove_device_obj_t<V&&>, U&&>>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Vector V, Scalar U>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto divide(U&& x, V&& v) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Div, U&&, remove_device_obj_t<V&&>>>(std::forward<U>(x), std::forward<V>(v));
    }

    template<Vector V1, Vector V2>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto divide(V1&& v1, V2&& v2) noexcept requires(DeviceObj<V1> && DeviceObj<V2>) {
        return device_obj<VectorExpr<ExprID::Div, remove_device_obj_t<V1&&>, remove_device_obj_t<V2&&>>>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
