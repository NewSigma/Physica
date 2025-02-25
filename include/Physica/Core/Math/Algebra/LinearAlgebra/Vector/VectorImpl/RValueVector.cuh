/*
 * Copyright 2022-2025 Weibo He.
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

#include "Physica/Core/Utils/CUDA/device_obj.cuh"
#include "RValueVector.h"

namespace Physica {
    template<class Derived>
    class device_obj<RValueVector<Derived>> : public CRTPBase<device_obj<RValueVector<Derived>>> {
        using This = device_obj<RValueVector<Derived>>;
        using Base = CRTPBase<This>;
        using TraitsType = Traits<device_obj<Derived>>;
    public:
        using ScalarType = TraitsType::ScalarType;
        using ValueType = ScalarType::ValueType;
        constexpr static size_t SizeAtCompile = TraitsType::SizeAtCompile;
        using PacketType = device_obj<BestPacket<ScalarType, SizeAtCompile>>::Type;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isComplex = ScalarType::isComplex;
        constexpr static size_t MaxThreadPerBlock = 256;
    protected:
        using T = ScalarType;
        using Tr = ScalarType::RealType;
    public:
        ~device_obj() = default;
        /* Operations */
        template<Vector V>
        __host__ __device__ void assign(device_obj<V>& target) const;

        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return Base::getDerived().calc(index); }
        template<Packet Pack>
        [[nodiscard]] __device__ inline Pack packet(size_t index) const;
        template<Packet Pack>
        [[nodiscard]] __device__ inline Pack packetPartial(size_t index, size_t count) const;

        [[nodiscard]] __host__ __device__ inline auto transpose() const noexcept;

        [[nodiscard]] __device__ inline Tr norm() const;
        [[nodiscard]] __device__ inline Tr squaredNorm() const;
        [[nodiscard]] __device__ T max() const;
        [[nodiscard]] __device__ T min() const;
        [[nodiscard]] __device__ T sum() const;
        template<Vector V>
        [[nodiscard]] __device__ inline auto crossProduct(const device_obj<V>& v) const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Base::getDerived().getLength(); }
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        template<Vector V>
        __device__ void assignToImpl(device_obj<V>& target) const;
    };
}

namespace Physica {
    template<class T>
    class Traits<device_obj<RValueVector<T>>> {
    public:
        using Derived = device_obj<T>;
    };
}

#include "RValueVectorImpl/RValueVectorImpl.cuh"
#include "RValueVectorImpl/CrossProduct.cuh"
#include "RValueVectorImpl/InnerDot.cuh"
#include "VectorExpr.cuh"

