/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <iostream>
#include <cmath>
#include <Physica/Core/Utils/DebugUtil.h>

namespace Physica {
    namespace Test {
        void checkGPU() {
            std::cout << "Checking for CUDA devices.\n";
            int deviceCount = 0;
            cudaGetDeviceCount(&deviceCount);
            if (deviceCount == 0) {
                std::cout << "No CUDA devices detected." << std::endl;
                return;
            }
            std::cout << "Detected " << deviceCount << " CUDA Capable device(s)\n";

            int driverVersion = 0, runtimeVersion = 0;
            cudaDeviceProp deviceProp{};
            for(int index = 0; index < deviceCount; ++index) {
                checkCudaError(cudaSetDevice(index));
                checkCudaError(cudaGetDeviceProperties(&deviceProp, index));
                cudaDriverGetVersion(&driverVersion);
                cudaRuntimeGetVersion(&runtimeVersion);

                std::cout << "Device " << index << ": \"" << deviceProp.name << "\"\n";

                std::cout << "  CUDA Driver Version / Runtime Version          "
                          << driverVersion / 1000 << '.' << (driverVersion % 100) / 10 << '/'
                          << runtimeVersion / 1000 << '.' << (runtimeVersion % 100) / 10 << '\n';
                std::cout << "  CUDA Capability Major/Minor version number:    "
                          << deviceProp.major << '.' << deviceProp.minor << '\n';
                std::cout << "  Total amount of global memory:                 "
                          << (float)deviceProp.totalGlobalMem / pow(1024.0, 3) << "GBytes ("
                          << (unsigned long long)deviceProp.totalGlobalMem << " bytes)\n";
                std::cout << "  GPU Clock rate:                                "
                          << (float)deviceProp.clockRate * 1e-3f << " MHz ("
                          << (float)deviceProp.clockRate * 1e-6f << " GHz)\n";
                std::cout << "  Memory Clock rate:                             "
                          << (float)deviceProp.memoryClockRate * 1e-3f << " Mhz\n";
                std::cout << "  Memory Bus Width:                              "
                          << deviceProp.memoryBusWidth << "-bit\n";

                if (deviceProp.l2CacheSize)
                    std::cout << "  L2 Cache Size:                                 " << deviceProp.l2CacheSize << " bytes\n";
                else
                    std::cout << "  L2 Cache Size:                                 None\n";

                std::cout << "  Max Texture Dimension Size (x,y,z)             "
                          << "1D=(" << deviceProp.maxTexture1D << "), "
                          << "2D=(" << deviceProp.maxTexture2D[0] << ',' << deviceProp.maxTexture2D[1] << "), "
                          << "3D=(" << deviceProp.maxTexture3D[0] << ',' << deviceProp.maxTexture3D[1] << ',' << deviceProp.maxTexture3D[2]  << ")\n";
                std::cout << "  Max Layered Texture Size (dim) x layers        "
                          << "1D=(" << deviceProp.maxTexture1DLayered[0] << ") x " << deviceProp.maxTexture1DLayered[1] << ", "
                          << "2D=(" << deviceProp.maxTexture2DLayered[0] << ", " << deviceProp.maxTexture2DLayered[1] << ") x " << deviceProp.maxTexture2DLayered[2] << '\n';
                std::cout << "  Total amount of constant memory:               " << deviceProp.totalConstMem << " bytes\n";
                std::cout << "  Total amount of shared memory per block:       " << deviceProp.sharedMemPerBlock << " bytes\n";
                std::cout << "  Total number of registers available per block: " << deviceProp.regsPerBlock << '\n';
                std::cout << "  Warp size:                                     " << deviceProp.warpSize << '\n';
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
                std::cout << "  Maximum memory pitch:                          " <<deviceProp.memPitch << " bytes\n" << std::endl;
            }
        }
    }
}
