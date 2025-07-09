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

#include "../VectorExpr.cuh"

namespace Physica {
    template<Vector V, Scalar U>
    class device_obj<VectorExpr<ExprType::Add, V, U>>
            : public device_obj<BinaryVectorExpr<ExprType::Add, V, U>> {
        using Base = device_obj<BinaryVectorExpr<ExprType::Add, V, U>>;
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
        void reverse(const Vector auto& grad) const noexcept requires(isReverseDiff) {
            const auto& g = grad.values();
            if constexpr (ReverseDiff<T>)
                Base::getLHS().reverse(g);
            if constexpr (ReverseDiff<U>)
                Base::getRHS().reverse(g.sum());
        }
    };

    template<Vector T1, Vector T2>
    class device_obj<VectorExpr<ExprType::Add, T1, T2>>
            : public device_obj<BinaryVectorExpr<ExprType::Add, T1, T2>> {
        using Base = device_obj<BinaryVectorExpr<ExprType::Add, T1, T2>>;
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
        void reverse(const Vector auto& grad) const noexcept requires(isReverseDiff);
        auto values() const noexcept { return Base::getLHS().values() + Base::getRHS().values(); }
    };

    template<Vector T1, Vector T2>
    void device_obj<VectorExpr<ExprType::Add, T1, T2>>::reverse(const Vector auto& grad) const noexcept requires(isReverseDiff) {
        const auto& g = grad.values();
        assert(g.getLength() == Base::getLength());
        if constexpr (ReverseDiff<T1>)
            Base::getLHS().reverse(g);
        if constexpr (ReverseDiff<T2>)
            Base::getRHS().reverse(g);
    }

    template<Vector T, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator+(T&& v, U&& x) noexcept requires(CUDA<T>) {
        return device_obj<VectorExpr<ExprType::Add, T&&, U&&>>(std::forward<T>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector T>
    [[nodiscard]] __host__ __device__ auto operator+(U&& x, T&& v) noexcept requires(CUDA<T>) {
        return std::forward<T>(v) + std::forward<U>(x);
    }

    template<Vector T1, Vector T2>
    [[nodiscard]] __host__ __device__ auto operator+(T1&& v1, T2&& v2) noexcept requires(CUDA<T1> && CUDA<T2>) {
        return device_obj<VectorExpr<ExprType::Add, T1&&, T2&&>>(std::forward<T1>(v1), std::forward<T2>(v2));
    }
}
