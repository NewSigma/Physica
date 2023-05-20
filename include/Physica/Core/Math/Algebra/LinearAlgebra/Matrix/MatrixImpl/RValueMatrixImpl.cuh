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

#include "Physica/Core/Parallel/StreamPool.cuh"

namespace Physica::Core {
    namespace Internal {
        template<class Derived, class OtherDerived>
        __global__ void assignTo_kernel(device_obj<RValueMatrix<Derived>> source, device_obj<LValueMatrix<OtherDerived>> target) {
            const size_t major = blockIdx.y;
            const size_t minor = blockIdx.x * blockDim.x + threadIdx.x;
            if (minor < source.getMaxMinor())
                target.refFromMajorMinor(major, minor) = source.calcFromMajorMinor(major, minor);
        }
    }

    template<class Derived>
    template<class OtherDerived>
    void device_obj<RValueMatrix<Derived>>::assignTo(device_obj<LValueMatrix<OtherDerived>>& target) const {
        constexpr int elemPerThread = 64;
        int device;
        cudaGetDevice(&device);
        const int maxThreadsPerBlock = Utils::DeviceProp::getInstance().getProperty(device).maxThreadsPerBlock;
        const size_t major = getMaxMajor();
        const size_t minor = getMaxMinor();
        const size_t numThread = minor > maxThreadPerBlock ? maxThreadPerBlock : minor;
        const size_t numBlockX = (minor + maxThreadsPerBlock) / maxThreadsPerBlock;
        Internal::assignTo_kernel<<<{numBlockX, major}, numThread, 0, StreamPool::getStream()>>>(Base::getDerived(), target.getDerived());
    }

    template<class Derived>
    typename device_obj<RValueMatrix<Derived>>::ScalarType
    __device__ device_obj<RValueMatrix<Derived>>::calcFromMajorMinor(size_t major, size_t minor) const {
        return calc(MatrixOption::rowFromMajorMinor<Derived>(major, minor), MatrixOption::columnFromMajorMinor<Derived>(major, minor));
    }
}
