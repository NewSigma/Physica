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

namespace Physica::Core {
    namespace Internal {
        template<class Derived, class OtherDerived>
        __global__ void assignTo_kernel(
                Physica::PlainStruct<Derived> source,
                Physica::PlainStruct<OtherDerived> target) {
            using namespace Physica::Core;
            using HostDerived = typename Derived::host_obj;
            using HostOtherDerived = typename OtherDerived::host_obj;
            static_assert(std::is_base_of<RValueVector<HostDerived>, HostDerived>::value, "[Error]: Invalid source vector type");
            static_assert(std::is_base_of<LValueVector<HostOtherDerived>, HostOtherDerived>::value, "[Error]: Invalid target vector type");
            const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = source.getDerived().getLength();
            if (index < length)
                target.getDerived()[index] = source.getDerived().calc(index);
        }
    }

    template<class Derived>
    template<class OtherDerived>
    void device_obj<RValueVector<Derived>>::assignTo(device_obj<LValueVector<OtherDerived>>& target) const {
        using namespace Physica;
        int device;
        cudaGetDevice(&device);
        const int maxThreadsPerBlock = Utils::DeviceProp::getInstance().getProperty(device).maxThreadsPerBlock;
        const int numBlock = (getLength() + maxThreadsPerBlock) / maxThreadsPerBlock;
        const int numThread = getLength() >= maxThreadsPerBlock ? maxThreadsPerBlock : getLength();
        Internal::assignTo_kernel<<<numBlock, numThread>>>(asStruct(Base::getDerived()), asStruct(target.getDerived()));
    }
}
