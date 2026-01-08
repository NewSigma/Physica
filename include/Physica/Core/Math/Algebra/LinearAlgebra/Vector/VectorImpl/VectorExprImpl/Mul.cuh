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
    class device_obj<VectorExpr<ExprID::Mul, V, U>>
            : public device_obj<BinaryVectorExpr<ExprID::Mul, V, U>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Mul, V, U>>;
    public:
        using typename Base::ScalarType;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            if constexpr (isReverseDiff)
                return calc_value(index);
            else
                return Base::getLHS().calc(index) * Base::getRHS();
        }

        [[nodiscard]] __device__ Tv calc_value(size_t index) const {
            return Base::getLHS().calc_value(index) * Base::getRHS().value();
        }

        template<Packet Pack>
        [[nodiscard]] __device__ Pack packet(size_t index) const {
            return Base::getLHS().template packet<Pack>(index) * Pack(Base::getRHS());
        }

        template<Packet Pack>
        [[nodiscard]] __device__ Pack packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<Pack>(index, count) * Pack(Base::getRHS());
        }
    };

    template<Vector V1, Vector V2>
    class device_obj<VectorExpr<ExprID::Mul, V1, V2>>
            : public device_obj<BinaryVectorExpr<ExprID::Mul, V1, V2>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Mul, V1, V2>>;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t index) const {
            return Base::getLHS().calc(index) * Base::getRHS().calc(index);
        }

        template<Packet Pack>
        [[nodiscard]] __device__ Pack packet(size_t index) const {
            return Base::getLHS().template packet<Pack>(index)
                 * Base::getRHS().template packet<Pack>(index);
        }

        template<Packet Pack>
        [[nodiscard]] __device__ Pack packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<Pack>(index, count)
                 * Base::getRHS().template packetPartial<Pack>(index, count);
        }
    };

    template<Vector V, Scalar U>
    [[nodiscard]] __host__ __device__ auto operator*(V&& v, U&& x) noexcept requires(CUDA<V>) {
        return device_obj<VectorExpr<ExprID::Mul, V&&, U&&>>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard]] __host__ __device__ auto operator*(U&& x, V&& v) noexcept requires(CUDA<V>) {
        return std::forward<V>(v) * std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] __host__ __device__ auto hadamard(V1&& v1, V2&& v2) noexcept requires(CUDA<V1> && CUDA<V2>) {
        return device_obj<VectorExpr<ExprID::Mul, V1&&, V2&&>>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
