/*
 * Copyright 2022-2024 Weibo He.
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
#include "Physica/Core/Utils/CUDA/DeviceProp.cuh"
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
    private:
        using RealType = ScalarType::RealType;
    public:
        ~device_obj() = default;
        /* Operations */
        template<Vector V>
        __host__ __device__ void assign(device_obj<V>& target) const;
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return Base::getDerived().template calc<Owner>(index); }
        template<class AnyPacket, Side Owner = GetSide()>
        [[nodiscard]] __device__ inline AnyPacket packet(size_t index) const;
        template<class AnyPacket, Side Owner = GetSide()>
        [[nodiscard]] __device__ inline AnyPacket packetPartial(size_t index, size_t count) const;

        [[nodiscard]] __host__ __device__ inline auto transpose() const noexcept;
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Base::getDerived().template getLength<Owner>(); }
        [[nodiscard]] __device__ inline RealType norm() const;
        [[nodiscard]] __device__ inline RealType squaredNorm() const;
        [[nodiscard]] __device__ ScalarType max() const;
        [[nodiscard]] __device__ ScalarType min() const;
        [[nodiscard]] __device__ ScalarType sum() const;
        template<Vector V>
        [[nodiscard]] __device__ inline device_obj<CrossProduct<Derived, V>> crossProduct(const device_obj<V>& v) const noexcept;
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

