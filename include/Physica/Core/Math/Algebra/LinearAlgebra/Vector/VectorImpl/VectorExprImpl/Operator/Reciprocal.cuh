/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.cuh"

namespace Physica {
    template<Vector V>
    class device_obj<VectorExpr<ExprID::Reciprocal, V>>
            : public device_obj<UnaryVectorExpr<ExprID::Reciprocal, V>> {
        using Base = device_obj<UnaryVectorExpr<ExprID::Reciprocal, V>>;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Getters */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;

        template<int Size>
        [[nodiscard]] __device__ device_obj<SIMD<T, Size>> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] __device__ device_obj<SIMD<T, Size>> packet(size_t index, size_t count) const noexcept;
    };

    template<Vector V>
    __device__ auto device_obj<VectorExpr<ExprID::Reciprocal, V>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        return reciprocal(Base::getExpr().calc(index, block));
    }

    template<Vector V>
    template<int Size>
    __device__ auto device_obj<VectorExpr<ExprID::Reciprocal, V>>::packet(size_t index) const noexcept -> device_obj<SIMD<T, Size>> {
        return reciprocal(Base::getExpr().template packet<Size>(index));
    }

    template<Vector V>
    template<int Size>
    __device__ auto device_obj<VectorExpr<ExprID::Reciprocal, V>>::packet(size_t index, size_t count) const noexcept -> device_obj<SIMD<T, Size>> {
        return reciprocal(Base::getExpr().template packet<Size>(index, count));
    }

    template<Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto reciprocal(V&& v) noexcept requires(DeviceObj<V> && Vector<V>) {
        return device_obj<VectorExpr<ExprID::Reciprocal, remove_device_obj_t<V&&>>>(std::forward<V>(v));
    }
}
