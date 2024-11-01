/*
 * Copyright 2019-2024 Weibo He.
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

#include <climits>
#include <cstddef>
#include <Physica/Config.h>

//Avoid "unused parameter" warnings and note that the parameter is used in asm.
#define USE_IN_ASM(x) (void)x
/**
 * Improve: Platform dependent, may be wrong on some spatial platforms, add a test before compile.
 * Note: Use of double_extract may cause several warnings in valgrind.
 */
#if PhysicaEndian == PhysicaBigEndian
    union float_extract {
        float value;
        struct {
            unsigned int sign : 1;
            unsigned int exp : 8;
            unsigned int high : 7;
            unsigned int low : 16;
        };
    };

    union double_extract {
        double value;
        struct {
            unsigned int sign : 1;
            unsigned int exp : 11;
            unsigned int high : 20;
            unsigned int low : 32;
        };
    };
#elif PhysicaEndian == PhysicaLittleEndian
    union float_extract {
        float value;
        struct {
            unsigned int low : 16;
            unsigned int high : 7;
            unsigned int exp : 8;
            unsigned int sign : 1;
        };
    };

    union double_extract {
        double value;
        struct {
            unsigned int low : 32;
            unsigned int high : 20;
            unsigned int exp : 11;
            unsigned int sign : 1;
        };
    };
#endif

#ifdef PHYSICA_MKL
    #include <mkl_types.h>
#endif

#ifndef PHYSICA_CUDA
    #define __host__
    #define __device__
#else
    #include <cuda_runtime_api.h>
#endif

#ifdef __GNUC__
    #define PHYSICA_API __attribute__((visibility("default")))
    #ifndef __forceinline
        #define __forceinline inline __attribute__((always_inline))
    #endif
#else
    #define PHYSICA_API __declspec(dllexport)
    #define __restrict__ __restrict

    #include <BaseTsd.h>
    using ssize_t = SSIZE_T;
#endif

namespace Physica {
    template<class T> class Traits;
    /**
     * Forward declaration for friend-class-based tests
     */
    class Test;

    constexpr static unsigned int PhysicaWordSize = sizeof(void*) * CHAR_BIT;
    constexpr size_t Dynamic = 0;

    enum class Side {
        Host,
        Device
    };

    __host__ __device__ constexpr inline static Side GetSide() {
    #ifdef __CUDA_ARCH__
        return Side::Device;
    #else
        return Side::Host;
    #endif
    }

    __host__ __device__ constexpr inline static bool IsHost() {
        return GetSide() == Side::Host;
    }

    __host__ __device__ constexpr inline static bool IsDevice() {
        return GetSide() == Side::Device;
    }

    constexpr inline static bool IsMSVC() {
    #ifdef _MSC_VER
        return true;
    #else
        return false;
    #endif
    }

    constexpr inline static bool UseASM() {
    #ifdef PHYSICA_ASM
        return true;
    #else
        return false;
    #endif
    }

    constexpr inline static bool HasHDF5() {
    #ifdef PHYSICA_HDF5
        return true;
    #else
        return false;
    #endif
    }

    constexpr inline static bool HasMKL() {
    #ifdef PHYSICA_MKL
        return true;
    #else
        return false;
    #endif
    }

    constexpr inline static bool HasMPI() {
    #ifdef PHYSICA_MPI
        return true;
    #else
        return false;
    #endif
    }

    constexpr inline static bool HasCUDA() {
    #ifdef PHYSICA_CUDA
        return true;
    #else
        return false;
    #endif
    }
}
