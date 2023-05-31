/*
 * Copyright 2022 WeiBo He.
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
        __device__ inline void assignToImpl(const Derived& source, OtherDerived& target) {
            const unsigned int delta = gridDim.x * blockDim.x;
            const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = source.getLength();
            for (unsigned int shift = 0; shift < length; shift += delta) {
                const unsigned int index = id + shift;
                if (index < length)
                    target[index] = source.calc(index);
            }
        }

        template<class Derived, class OtherDerived>
        __global__ void assignTo_kernel(
                Physica::PlainStruct<Derived> source,
                Physica::PlainStruct<OtherDerived> target) {
            using namespace Physica::Core;
            using HostDerived = typename Derived::host_obj;
            using HostOtherDerived = typename OtherDerived::host_obj;
            static_assert(std::is_base_of<RValueVector<HostDerived>, HostDerived>::value, "[Error]: Invalid source vector type");
            static_assert(std::is_base_of<LValueVector<HostOtherDerived>, HostOtherDerived>::value, "[Error]: Invalid target vector type");
            assignToImpl(source.getDerived(), target.getDerived());
        }
    }

    template<class Derived>
    template<class OtherDerived>
    __host__ __device__
    void device_obj<RValueVector<Derived>>::assignTo(device_obj<LValueVector<OtherDerived>>& target) const {
        [[maybe_unused]] const auto _ = Internal::assignTo_kernel<device_obj<Derived>, device_obj<OtherDerived>>;
    #ifndef  __CUDA_ARCH__
        using namespace Physica;
        constexpr unsigned int WarpSize = Utils::DeviceProp::WarpSize;
        const int numBlock = (getLength() + WarpSize - 1) / WarpSize;
        const int numThread = WarpSize;
        Internal::assignTo_kernel<<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(Base::getDerived()), asStruct(target.getDerived()));
    #else
        Internal::assignToImpl(Base::getDerived(), target.getDerived());
    #endif
    }
}
