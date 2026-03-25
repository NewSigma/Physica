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
    template<Vector V, Scalar U>
    class device_obj<VectorExpr<ExprID::Add, V, U>>
            : public device_obj<BinaryVectorExpr<ExprID::Add, V, U>> {
        using host_obj = VectorExpr<ExprID::Add, V, U>;
        using Base = device_obj<BinaryVectorExpr<ExprID::Add, V, U>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t index) const;
        [[nodiscard]] __device__ Tv calc_value(size_t index) const;

        template<int Size>
        [[nodiscard]] __device__ SIMD<T, Size> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] __device__ SIMD<T, Size> packet(size_t index, size_t count) const noexcept;

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector V, Scalar U>
    [[nodiscard]] __device__ auto device_obj<VectorExpr<ExprID::Add, V, U>>::calc(size_t index) const -> T {
        if constexpr (host_obj::lowerToFMA()) {
            T a = getLHS().getLHS().calc(index);
            T b;
            if constexpr (Scalar<decltype(getLHS().getRHS())>)
                b = getLHS().getRHS();
            else
                b = getLHS().getRHS().calc(index);
            T c = getRHS();
            return fma(a, b, c);
        }
        else if constexpr (isReverseDiff)
            return calc_value(index);
        else
            return Base::getLHS().calc(index) + Base::getRHS();
    }

    template<Vector V, Scalar U>
    [[nodiscard]] __device__ auto device_obj<VectorExpr<ExprID::Add, V, U>>::calc_value(size_t index) const -> Tv {
        return Base::getLHS().calc_value(index) + Base::getRHS().value();
    }

    template<Vector V, Scalar U>
    template<int Size>
    [[nodiscard]] __device__ auto device_obj<VectorExpr<ExprID::Add, V, U>>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        return Base::getLHS().template packet<Size>(index) + SIMD<T, Size>(Base::getRHS());
    }

    template<Vector V, Scalar U>
    template<int Size>
    [[nodiscard]] __device__ auto device_obj<VectorExpr<ExprID::Add, V, U>>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size>{
        return Base::getLHS().template packet<Size>(index, count) + SIMD<T, Size>(Base::getRHS());
    }

    template<Vector V, Scalar U>
    void device_obj<VectorExpr<ExprID::Add, V, U>>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& g = grad.values();
        if constexpr (ReverseDiff<V>)
            Base::getLHS().reverse(g);
        if constexpr (ReverseDiff<U>)
            Base::getRHS().reverse(g.sum());
    }

    template<Vector V1, Vector V2>
    class device_obj<VectorExpr<ExprID::Add, V1, V2>>
            : public device_obj<BinaryVectorExpr<ExprID::Add, V1, V2>> {
        using host_obj = VectorExpr<ExprID::Add, V1, V2>;
        using Base = device_obj<BinaryVectorExpr<ExprID::Add, V1, V2>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t index) const;
        [[nodiscard]] __device__ Tv calc_value(size_t index) const;

        template<int Size>
        [[nodiscard]] __device__ SIMD<T, Size> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] __device__ SIMD<T, Size> packet(size_t index, size_t count) const noexcept;

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;
        auto values() const noexcept { return Base::getLHS().values() + Base::getRHS().values(); }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    };

    template<Vector V1, Vector V2>
    __device__ auto device_obj<VectorExpr<ExprID::Add, V1, V2>>::calc(size_t index) const -> T {
        if constexpr (host_obj::lowerToFMA()) {
            T a = getLHS().getLHS().calc(index);
            T b;
            if constexpr (Scalar<decltype(getLHS().getRHS())>)
                b = getLHS().getRHS();
            else
                b = getLHS().getRHS().calc(index);
            T c = getRHS().calc(index);
            return fma(a, b, c);
        }
        else if constexpr (isReverseDiff)
            return calc_value(index);
        else
            return Base::getLHS().calc(index) + Base::getRHS().calc(index);
    }

    template<Vector V1, Vector V2>
    __device__ auto device_obj<VectorExpr<ExprID::Add, V1, V2>>::calc_value(size_t index) const -> Tv {
        return Base::getLHS().calc_value(index) + Base::getRHS().calc_value(index);
    }

    template<Vector V1, Vector V2>
    template<int Size>
    __device__ auto device_obj<VectorExpr<ExprID::Add, V1, V2>>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        return Base::getLHS().template packet<Size>(index)
             + Base::getRHS().template packet<Size>(index);
    }

    template<Vector V1, Vector V2>
    template<int Size>
    __device__ auto device_obj<VectorExpr<ExprID::Add, V1, V2>>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size> {
        return Base::getLHS().template packet<Size>(index, count)
             + Base::getRHS().template packet<Size>(index, count);
    }

    template<Vector V1, Vector V2>
    void device_obj<VectorExpr<ExprID::Add, V1, V2>>::reverse(const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        const auto& g = grad.values();
        assert(g.getLength() == Base::getLength());
        if constexpr (ReverseDiff<V1>)
            Base::getLHS().reverse(g);
        if constexpr (ReverseDiff<V2>)
            Base::getRHS().reverse(g);
    }

    template<Vector V, Scalar U>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator+(V&& v, U&& x) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Add, remove_device_obj_t<V&&>, U&&>>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator+(U&& x, V&& v) noexcept requires(DeviceObj<V>) {
        return std::forward<V>(v) + std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator+(V1&& v1, V2&& v2) noexcept requires(DeviceObj<V1> && DeviceObj<V2>) {
        if constexpr (!canonicalized(v1, v2))
            return std::forward<V2>(v2) + std::forward<V1>(v1);
        else
            return device_obj<VectorExpr<ExprID::Add, remove_device_obj_t<V1&&>, remove_device_obj_t<V2&&>>>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
