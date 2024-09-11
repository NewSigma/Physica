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

#include <Physica/Utils/CUDA/device_obj.cuh>
#include <Physica/Utils/CUDA/DeviceProp.cuh>
#include "RValueVector.h"

namespace Physica::Core {
    template<class T>
    struct is_vector<device_obj<T>> : public is_vector<T> {};

    template<class Derived>
    class device_obj<RValueVector<Derived>> : public Utils::CRTPBase<device_obj<Derived>> {
        using Base = Utils::CRTPBase<device_obj<Derived>>;
        using TraitsType = Traits<device_obj<Derived>>;
    public:
        using ScalarType = typename TraitsType::ScalarType;
        constexpr static size_t SizeAtCompile = TraitsType::SizeAtCompile;
        using PacketType = typename device_obj<BestPacket<ScalarType, SizeAtCompile>>::Type;
        constexpr static bool isComplex = ScalarType::isComplex;
        constexpr static size_t MaxThreadPerBlock = 256;
    private:
        using RealType = typename ScalarType::RealType;
    public:
        ~device_obj() = default;
        /* Operations */
        template<class OtherDerived>
        __host__ __device__ void assignTo(device_obj<LValueVector<OtherDerived>>& target) const;
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return Base::getDerived().template calc<Owner>(index); }
        [[nodiscard]] __host__ __device__ inline device_obj<TransposeVector<Derived>> transpose() const noexcept;
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Base::getDerived().template getLength<Owner>(); }
        [[nodiscard]] __device__ inline RealType norm() const;
        [[nodiscard]] __device__ inline RealType squaredNorm() const;
        [[nodiscard]] __device__ ScalarType max() const;
        [[nodiscard]] __device__ ScalarType min() const;
        [[nodiscard]] __device__ ScalarType sum() const;
        template<class OtherDerived>
        [[nodiscard]] __device__ inline device_obj<CrossProduct<Derived, OtherDerived>> crossProduct(const device_obj<RValueVector<OtherDerived>>& v) const noexcept;
    protected:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = default;
        device_obj& operator=(device_obj&&) noexcept = default;
        /* Operations */
        template<class OtherDerived>
        __device__ void assignToImpl(device_obj<LValueVector<OtherDerived>>& target_) const;
    };
}

#include "RValueVectorImpl/RValueVectorImpl.cuh"
#include "VectorExpr.cuh"
#include "InnerDot.cuh"

