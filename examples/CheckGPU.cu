/*
 * Copyright 2019 WeiBo He.
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
#include <cmath>
#include <Physica/Utils/CUDA/DebugUtil.cuh>

int main(int argc, char** argv) {
    std::cout << "Checking for CUDA devices.\n";
    /* Version */ {
        int driverVersion = 0, runtimeVersion = 0;
        cudaDriverGetVersion(&driverVersion);
        cudaRuntimeGetVersion(&runtimeVersion);
        std::cout << "  CUDA Driver Version / Runtime Version          "
                  << driverVersion / 1000 << '.' << (driverVersion % 100) / 10 << '/'
                  << runtimeVersion / 1000 << '.' << (runtimeVersion % 100) / 10 << '\n';
    }
    //Get device count.
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0) {
        std::cout << "No CUDA devices detected." << std::endl;
        return 1;
    }
    std::cout << "Detected " << deviceCount << " CUDA Capable device(s)\n";
    //Print properties for all devices.
    cudaDeviceProp deviceProp{};
    for(int index = 0; index < deviceCount; ++index) {
        checkCudaError(cudaGetDeviceProperties(&deviceProp, index));

        std::cout << "Device " << index << ": \"" << deviceProp.name << "\"\n";

        std::cout << "  CUDA Capability Major/Minor version number:    "
                  << deviceProp.major << '.' << deviceProp.minor << '\n';
        /* Memory amount */
        std::cout << "  Memory amount:\n";
        std::cout << "  Total amount of global memory:                 "
                  << (float)deviceProp.totalGlobalMem / pow(1024.0, 3) << "GBytes ("
                  << (unsigned long long)deviceProp.totalGlobalMem << " bytes)\n";
        std::cout << "  Shared memory per multiprocessor:              " << deviceProp.sharedMemPerMultiprocessor << " bytes\n";
        std::cout << "  Shared memory per block:                       " << deviceProp.sharedMemPerBlock << " bytes\n";
        std::cout << "  Total amount of constant memory:               " << deviceProp.totalConstMem << " bytes\n";
        std::cout << "  Maximum memory pitch:                          " <<deviceProp.memPitch << " bytes\n";
        /* Registers */
        std::cout << "  Registers amount:\n";
        std::cout << "  Registers available per multiprocessor:        " << deviceProp.regsPerMultiprocessor << '\n';
        std::cout << "  Registers available per block:                 " << deviceProp.regsPerBlock << '\n';
        std::cout << "  Warp size:                                     " << deviceProp.warpSize << '\n';
        /* Dimension */
        std::cout << "  Dimension:\n";
        std::cout << "  Maximum number of threads per multiprocessor:  " << deviceProp.maxThreadsPerMultiProcessor << '\n';
        std::cout << "  Maximum number of threads per block:           " << deviceProp.maxThreadsPerBlock << '\n';
        std::cout << "  Maximum sizes of each dimension of a block:    "
                  << deviceProp.maxThreadsDim[0] << " x "
                  << deviceProp.maxThreadsDim[1] << " x "
                  << deviceProp.maxThreadsDim[2] << '\n';
        std::cout << "  Maximum sizes of each dimension of a grid:     "
                  << deviceProp.maxGridSize[0] << " x "
                  << deviceProp.maxGridSize[1] << " x "
                  << deviceProp.maxGridSize[2] << '\n';
        std::cout << "  Max Texture Dimension Size (x,y,z)             "
                  << "1D=(" << deviceProp.maxTexture1D << "), "
                  << "2D=(" << deviceProp.maxTexture2D[0] << ',' << deviceProp.maxTexture2D[1] << "), "
                  << "3D=(" << deviceProp.maxTexture3D[0] << ',' << deviceProp.maxTexture3D[1] << ',' << deviceProp.maxTexture3D[2]  << ")\n";
        std::cout << "  Max Layered Texture Size (dim) x layers        "
                  << "1D=(" << deviceProp.maxTexture1DLayered[0] << ") x " << deviceProp.maxTexture1DLayered[1] << ", "
                  << "2D=(" << deviceProp.maxTexture2DLayered[0] << ", " << deviceProp.maxTexture2DLayered[1] << ") x " << deviceProp.maxTexture2DLayered[2] << '\n';
        /* Clock rates */
        std::cout << "  Clock rates:\n";
        std::cout << "  GPU Clock rate:                                "
                  << (float)deviceProp.clockRate * 1e-3f << " MHz ("
                  << (float)deviceProp.clockRate * 1e-6f << " GHz)\n";
        std::cout << "  Memory Clock rate:                             "
                  << (float)deviceProp.memoryClockRate * 1e-3f << " Mhz\n";
        std::cout << "  Memory Bus Width:                              "
                  << deviceProp.memoryBusWidth << "-bit\n";
        /* Caches */
        std::cout << "  Caches:\n";
        if (deviceProp.l2CacheSize)
            std::cout << "  L2 Cache Size:                                 " << deviceProp.l2CacheSize << " bytes\n";
        else
            std::cout << "  L2 Cache Size:                                 Null\n";
        std::cout << std::endl;
    }
    return 0;
}
