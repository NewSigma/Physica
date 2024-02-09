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

#include "Physica/Utils/CUDA/PlainStruct.h"
#include "Physica/Core/Parallel/StreamPool.cuh"

namespace Physica::Core {
    namespace Internal {
        template<class Derived, class OtherDerived>
        __global__ void RValueVector_assignToKernel(
                Physica::PlainStruct<const Derived> source,
                Physica::PlainStruct<OtherDerived> target) {
            static_assert(is_vector<Derived>::value, "[Error]: Invalid source vector type");
            static_assert(is_vector<OtherDerived>::value, "[Error]: Invalid target vector type");
            using namespace Physica::Core;
            using ScalarType = typename OtherDerived::ScalarType;

            const unsigned int delta = gridDim.x * blockDim.x;
            const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = source.getDerived().getLength();
            for (unsigned int shift = 0; shift < length; shift += delta) {
                const unsigned int index = id + shift;
                if (index < length)
                    target.getDerived()[index] = ScalarType(source.getDerived().calc(index));
            }
        }
    }

    template<class Derived>
    template<class OtherDerived>
    __host__ __device__
    void device_obj<RValueVector<Derived>>::assignTo(device_obj<LValueVector<OtherDerived>>& target) const {
        constexpr size_t OtherSize = Internal::Traits<OtherDerived>::SizeAtCompile;
        static_assert(SizeAtCompile == Dynamic || OtherSize == Dynamic || SizeAtCompile == OtherSize,
                "[Error]: Size mismatch between two vector");
        [[maybe_unused]] const auto _ = Internal::RValueVector_assignToKernel<device_obj<Derived>, device_obj<OtherDerived>>;
    #ifndef  __CUDA_ARCH__
        using namespace Physica;
        constexpr unsigned int WarpSize = Utils::DeviceProp::WarpSize;
        const int numBlock = (getLength() + WarpSize - 1) / WarpSize;
        const int numThread = WarpSize;
        Internal::RValueVector_assignToKernel<<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(Base::getDerived()), asStruct(target.getDerived()));
    #else
        assignToImpl(target.getDerived());
    #endif
    }

    template<class Derived>
    template<class OtherDerived>
    __device__ inline void device_obj<RValueVector<Derived>>::assignToImpl(OtherDerived& target) const {
        using OtherScalar = typename OtherDerived::ScalarType;
        for (size_t i = 0; i < getLength(); ++i)
            target[i] = OtherScalar(calc(i));
    }

    template<class Derived>
    __host__ __device__ inline device_obj<TransposeVector<Derived>> device_obj<RValueVector<Derived>>::transpose() const noexcept {
        return device_obj<TransposeVector<Derived>>(*this);
    }

    template<class Derived>
    __device__ inline typename device_obj<RValueVector<Derived>>::RealType device_obj<RValueVector<Derived>>::norm() const {
        return sqrt(Base::getDerived().squaredNorm());
    }

    template<class Derived>
    __device__ inline typename device_obj<RValueVector<Derived>>::RealType device_obj<RValueVector<Derived>>::squaredNorm() const {
        auto result = RealType(0);
        for(size_t i = 0; i < getLength(); ++i)
            result += calc(i).squaredNorm();
        return result;
    }

    template<class Derived>
    __device__ typename device_obj<RValueVector<Derived>>::ScalarType device_obj<RValueVector<Derived>>::max() const {
        assert(getLength() != 0);
        ScalarType result = calc(0);
        for(size_t i = 1; i < getLength(); ++i) {
            ScalarType temp = calc(i);
            if (result < temp)
                result = temp;
        }
        return result;
    }

    template<class Derived>
    __device__ typename device_obj<RValueVector<Derived>>::ScalarType device_obj<RValueVector<Derived>>::min() const {
        assert(getLength() != 0);
        ScalarType result = calc(0);
        for(size_t i = 1; i < getLength(); ++i) {
            ScalarType temp = calc(i);
            if (result > temp)
                result = temp;
        }
        return result;
    }

    template<class Derived>
    __device__ typename device_obj<RValueVector<Derived>>::ScalarType device_obj<RValueVector<Derived>>::sum() const {
        assert(getLength() != 0);
        auto result = ScalarType(0);
        for(size_t i = 0; i < getLength(); ++i)
            result += calc(i);
        return result;
    }

    template<class Derived>
    template<class OtherDerived>
    __device__ inline device_obj<CrossProduct<Derived, OtherDerived>>
    device_obj<RValueVector<Derived>>::crossProduct(const device_obj<RValueVector<OtherDerived>>& v) const noexcept {
        return device_obj<CrossProduct<Derived, OtherDerived>>(*this, v);
    }

    template<class VectorType1, class VectorType2>
    __device__ typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type
    operator*(const device_obj<RValueVector<VectorType1>>& v1, const device_obj<RValueVector<VectorType2>>& v2) {
        assert(v1.getLength() == v2.getLength());
        using ResultType = typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type;
        ResultType result = 0;
        for (size_t i = 0; i < v1.getLength(); ++i)
            result += v1.calc(i) * v2.calc(i);
        return result;
    }
}
