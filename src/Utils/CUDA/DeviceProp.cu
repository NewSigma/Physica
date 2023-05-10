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
#include <iostream>
#include "Physica/Utils/CUDA/DeviceProp.cuh"
#include "Physica/Core/Exception/CudaException.cuh"

namespace Physica::Utils {
    const DeviceProp& DeviceProp::getInstance() {
        static DeviceProp prop{};
        return prop;
    }

    DeviceProp::DeviceProp() {
        cudaError_t code;
        code = cudaDriverGetVersion(&driverVersion);
        if (code != cudaError_t::cudaSuccess)
            throw Core::CudaException(code);

        code = cudaRuntimeGetVersion(&runtimeVersion);
        if (code != cudaError_t::cudaSuccess)
            throw Core::CudaException(code);

        int deviceCount = 0;
        code = cudaGetDeviceCount(&deviceCount);
        if (code != cudaError_t::cudaSuccess)
            throw Core::CudaException(code);

        propList.resize(deviceCount);
        for (int i = 0; i < propList.getLength(); ++i) {
            code = cudaGetDeviceProperties(&propList[i], i);
            if (code != cudaError_t::cudaSuccess)
                throw Core::CudaException(code);
            if (getProperty(i).warpSize != WarpSize) {
                std::cerr << "[Error]: WarpSize is not standard\n";
                exit(EXIT_FAILURE);
            }
        }
    }
}
