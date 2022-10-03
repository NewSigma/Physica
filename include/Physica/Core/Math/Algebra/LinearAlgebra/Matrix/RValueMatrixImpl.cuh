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

namespace Physica::Core {
    namespace Internal {
        template<class Derived, class OtherDerived>
        __global__ void assignTo_kernel(device_obj<RValueMatrix<Derived>> source, device_obj<LValueMatrix<OtherDerived>> target) {
            //const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
            //if (index < getLength())
            //    target[index] = source.calc(index);
        }
    }

    template<class Derived>
    template<class OtherDerived>
    void device_obj<RValueMatrix<Derived>>::assignTo(device_obj<LValueMatrix<OtherDerived>>& target) const {
        //constexpr int elemPerThread = 64;
        //int device;
        //cudaGetDevice(&device);
        //const int maxThreadsPerBlock = Utils::DeviceProp::getInstance().getProperty(device).maxThreadsPerBlock;
        //const int numBlock = (getLength() + maxThreadsPerBlock) / maxThreadsPerBlock;
        //const int numThread = getLength() >= maxThreadPerBlock ? maxThreadPerBlock : getLength();
        //Internal::assignTo_kernel<<<numBlock, numThread>>>(Base::getDerived(), target.getDerived());
    }

    template<class Derived>
    typename device_obj<RValueMatrix<Derived>>::ScalarType
    device_obj<RValueMatrix<Derived>>::calcFromMajorMinor(size_t major, size_t minor) const {
        return calc(MatrixOption::rowFromMajorMinor<Derived>(major, minor), MatrixOption::columnFromMajorMinor<Derived>(major, minor));
    }
}
