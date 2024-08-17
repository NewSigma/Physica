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

#include "Physica/Utils/CUDA/PlainStruct.h"
#include "Physica/Core/Parallel/StreamPool.cuh"

namespace Physica::Core {
    namespace Internal {
        template<class Derived, class OtherDerived>
        __global__ void RValueVector_assignToKernel(
                Physica::PlainStruct<const Derived> source_,
                Physica::PlainStruct<OtherDerived> target_) {
            static_assert(is_vector<Derived>::value, "[Error]: Invalid source vector type");
            static_assert(is_vector<OtherDerived>::value, "[Error]: Invalid target vector type");
            using namespace Physica::Core;
            using ScalarType = typename OtherDerived::ScalarType;

            const auto& source = source_.getDerived();
            auto& target = target_.getDerived();
            const size_t length = source.getLength();
            const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
            if (index < length)
                target[index] = ScalarType(source.calc(index));
        }
    }

    template<class Derived>
    template<class OtherDerived>
    __host__ __device__
    void device_obj<RValueVector<Derived>>::assignTo(device_obj<LValueVector<OtherDerived>>& target) const {
        constexpr size_t OtherSize = Traits<OtherDerived>::SizeAtCompile;
        static_assert(SizeAtCompile == Dynamic || OtherSize == Dynamic || SizeAtCompile == OtherSize,
                "[Error]: Size mismatch between two vector");
        assert(target.getLength() == getLength() && "[Error]: Size mismatch between two vector");
        [[maybe_unused]] const auto kernel = Internal::RValueVector_assignToKernel<device_obj<Derived>, device_obj<OtherDerived>>;
    #ifndef  __CUDA_ARCH__
        using namespace Physica;
        const size_t length = getLength();
        const unsigned int numThread = std::min(length, MaxThreadPerBlock);
        const unsigned int numBlock = (length + numThread - 1) / numThread;
        kernel<<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(Base::getDerived()), asStruct(target.getDerived()));
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
}
