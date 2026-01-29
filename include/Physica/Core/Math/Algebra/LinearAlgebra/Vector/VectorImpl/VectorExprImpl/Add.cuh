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

#include "../VectorExpr.cuh"

namespace Physica {
    template<Vector V, Scalar U>
    class device_obj<VectorExpr<ExprID::Add, V, U>>
            : public device_obj<BinaryVectorExpr<ExprID::Add, V, U>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Add, V, U>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t index) const {
            if constexpr (isReverseDiff)
                return calc_value(index);
            else
                return Base::getLHS().calc(index) + Base::getRHS();
        }

        [[nodiscard]] __device__ Tv calc_value(size_t index) const {
            return Base::getLHS().calc_value(index) + Base::getRHS().value();
        }

        template<Packet Pack>
        [[nodiscard]] __device__ Pack packet(size_t index) const {
            return Base::getLHS().template packet<Pack>(index) + Pack(Base::getRHS());
        }

        template<Packet Pack>
        [[nodiscard]] __device__ Pack packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<Pack>(index, count) + Pack(Base::getRHS());
        }

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept {
            static_assert(isReverseDiff);
            const auto& g = grad.values();
            if constexpr (ReverseDiff<V>)
                Base::getLHS().reverse(g);
            if constexpr (ReverseDiff<U>)
                Base::getRHS().reverse(g.sum());
        }
    };

    template<Vector V1, Vector V2>
    class device_obj<VectorExpr<ExprID::Add, V1, V2>>
            : public device_obj<BinaryVectorExpr<ExprID::Add, V1, V2>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Add, V1, V2>>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t index) const {
            if constexpr (isReverseDiff)
                return calc_value(index);
            else
                return Base::getLHS().calc(index) + Base::getRHS().calc(index);
        }

        [[nodiscard]] __device__ Tv calc_value(size_t index) const {
            return Base::getLHS().calc_value(index) + Base::getRHS().calc_value(index);
        }

        template<Packet Pack>
        [[nodiscard]] __device__ Pack packet(size_t index) const {
            return Base::getLHS().template packet<Pack>(index)
                 + Base::getRHS().template packet<Pack>(index);
        }

        template<Packet Pack>
        [[nodiscard]] __device__ Pack packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<Pack>(index, count)
                 + Base::getRHS().template packetPartial<Pack>(index, count);
        }

        using Base::reverse;
        void reverse(const Vector auto& grad) const noexcept;
        auto values() const noexcept { return Base::getLHS().values() + Base::getRHS().values(); }
    };

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
    [[nodiscard]] __host__ __device__ auto operator+(V&& v, U&& x) noexcept requires(CUDA<V>) {
        return device_obj<VectorExpr<ExprID::Add, remove_device_obj_t<V&&>, U&&>>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard]] __host__ __device__ auto operator+(U&& x, V&& v) noexcept requires(CUDA<V>) {
        return std::forward<V>(v) + std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] __host__ __device__ auto operator+(V1&& v1, V2&& v2) noexcept requires(CUDA<V1> && CUDA<V2>) {
        return device_obj<VectorExpr<ExprID::Add, remove_device_obj_t<V1&&>, remove_device_obj_t<V2&&>>>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
