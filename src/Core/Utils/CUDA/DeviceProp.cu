/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <iostream>
#include "Physica/Core/Utils/CUDA/DeviceProp.cuh"
#include "Physica/Core/Exception/CUDA/CUDA.cuh"

namespace Physica {
    std::ostream& DeviceProp::printDeviceProp(unsigned int device, std::ostream& os) const {
        const auto& prop = getProperty(device);
        os << "Device " << device << ": \"" << prop.name << "\"\n";

        os << "  CUDA Capability Major/Minor version number:    "
                  << prop.major << '.' << prop.minor << '\n';
        /* Memory amount */
        os << "  Memory amount:\n";
        os << "  Total amount of global memory:                 "
                  << (float)prop.totalGlobalMem / pow(1024.0, 3) << " GBytes ("
                  << (unsigned long long)prop.totalGlobalMem << " bytes)\n";
        os << "  Shared memory per multiprocessor:              " << prop.sharedMemPerMultiprocessor << " bytes\n";
        os << "  Shared memory per block:                       " << prop.sharedMemPerBlock << " bytes\n";
        os << "  Total amount of constant memory:               " << prop.totalConstMem << " bytes\n";
        os << "  Maximum memory pitch:                          " << prop.memPitch << " bytes\n";
        /* Registers */
        os << "  Registers amount:\n";
        os << "  Registers available per multiprocessor:        " << prop.regsPerMultiprocessor << '\n';
        os << "  Registers available per block:                 " << prop.regsPerBlock << '\n';
        os << "  Warp size:                                     " << prop.warpSize << '\n';
        /* Dimension */
        os << "  Dimension:\n";
        os << "  Maximum number of threads per multiprocessor:  " << prop.maxThreadsPerMultiProcessor << '\n';
        os << "  Maximum number of threads per block:           " << prop.maxThreadsPerBlock << '\n';
        os << "  Maximum sizes of each dimension of a block:    "
                  << prop.maxThreadsDim[0] << " x "
                  << prop.maxThreadsDim[1] << " x "
                  << prop.maxThreadsDim[2] << '\n';
        os << "  Maximum sizes of each dimension of a grid:     "
                  << prop.maxGridSize[0] << " x "
                  << prop.maxGridSize[1] << " x "
                  << prop.maxGridSize[2] << '\n';
        os << "  Max Texture Dimension Size (x,y,z)             "
                  << "1D=(" << prop.maxTexture1D << "), "
                  << "2D=(" << prop.maxTexture2D[0] << ',' << prop.maxTexture2D[1] << "), "
                  << "3D=(" << prop.maxTexture3D[0] << ',' << prop.maxTexture3D[1] << ',' << prop.maxTexture3D[2]  << ")\n";
        os << "  Max Layered Texture Size (dim) x layers        "
                  << "1D=(" << prop.maxTexture1DLayered[0] << ") x " << prop.maxTexture1DLayered[1] << ", "
                  << "2D=(" << prop.maxTexture2DLayered[0] << ", " << prop.maxTexture2DLayered[1] << ") x " << prop.maxTexture2DLayered[2] << '\n';
        /* Clock rates */
        os << "  Clock rates:\n";
        os << "  GPU Clock rate:                                "
                  << (float)prop.clockRate * 1e-3f << " MHz ("
                  << (float)prop.clockRate * 1e-6f << " GHz)\n";
        os << "  Number of asynchronous engines:                " << prop.asyncEngineCount << '\n';
        os << "  Memory Clock rate:                             "
                  << (float)prop.memoryClockRate * 1e-3f << " Mhz\n";
        os << "  Memory Bus Width:                              "
                  << prop.memoryBusWidth << "-bit\n";
        /* Caches */
        os << "  Caches:\n";
        if (prop.l2CacheSize)
            os << "  L2 Cache Size:                                 " << prop.l2CacheSize << " bytes\n";
        else
            os << "  L2 Cache Size:                                 Null\n";
        os << std::endl;
        return os;
    }

    const DeviceProp& DeviceProp::getInstance() {
        static DeviceProp prop{};
        return prop;
    }

    DeviceProp::DeviceProp() {
        check(cudaDriverGetVersion(&driverVersion));
        check(cudaRuntimeGetVersion(&runtimeVersion));

        int deviceCount = 0;
        check(cudaGetDeviceCount(&deviceCount));
        propList.resize(deviceCount);
        for (size_t i = 0; i < propList.getLength(); ++i) {
            check(cudaGetDeviceProperties(&propList[i], i));
            if (getProperty(i).warpSize != Physica::CUDADevAttr::WarpSize) {
                std::cerr << "[Error]: WarpSize is not standard\n";
                exit(EXIT_FAILURE);
            }
        }
    }
}
