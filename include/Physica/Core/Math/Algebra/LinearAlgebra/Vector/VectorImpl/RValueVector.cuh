/*
 * Copyright 2022-2023 WeiBo He.
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

#include "RValueVector.h"
#include "Physica/Utils/CUDA/device_obj.cuh"
#include "Physica/Utils/CUDA/DeviceProp.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<RValueVector<Derived>> : public Utils::CRTPBase<device_obj<Derived>> {
        using Base = Utils::CRTPBase<device_obj<Derived>>;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
        constexpr static size_t SizeAtCompile = Internal::Traits<Derived>::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = Internal::Traits<Derived>::MaxSizeAtCompile;
        constexpr static bool isComplex = ScalarType::isComplex;
    private:
        using RealType = typename ScalarType::RealType;
    public:
        ~device_obj() = default;
        /* Operations */
        template<class OtherDerived>
        __host__ __device__ void assignTo(device_obj<LValueVector<OtherDerived>>& target) const;
        template<class OtherDerived>
        __device__ inline void assignToImpl(OtherDerived& target) const;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return Base::getDerived().calc(index); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Base::getDerived().getLength(); }
        [[nodiscard]] __device__ inline RealType norm() const;
        [[nodiscard]] __device__ inline RealType squaredNorm() const;
        [[nodiscard]] __device__ ScalarType max() const;
        [[nodiscard]] __device__ ScalarType min() const;
        [[nodiscard]] __device__ ScalarType sum() const;
    protected:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
    };

    template<class VectorType1, class VectorType2>
    __device__ typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type
    operator*(const device_obj<RValueVector<VectorType1>>& v1, const device_obj<RValueVector<VectorType2>>& v2);
}

#include "RValueVectorImpl.cuh"
#include "VectorExpression.cuh"

