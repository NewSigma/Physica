/*
 * Copyright 2019-2025 Weibo He.
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
#include "Physica/Config.h" // IWYU pragma: export
/**
 * \file Macro.h: Non-Proliferation of Macros
 */
#ifdef PHYSICA_MKL
    #include <mkl_types.h>
#else
    struct MKL_Complex8 {
        float real;
        float imag;
    };

    struct MKL_Complex16 {
        double real;
        double imag;
    };

    using VSLStreamStatePtr = void*;
#endif

#ifdef PHYSICA_CUDA
    #include <cuda_runtime.h>
#else
    #define __host__
    #define __device__
    struct curandGenerator_t {};
    using curandRngType_t = int;
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
    template<class T>
    class Traits;
    /**
     * Forward declaration for friend-class-based tests
     */
    class Test;

    [[maybe_unused]] constexpr static unsigned int PhysicaWordSize = sizeof(void*) * CHAR_BIT;
    [[maybe_unused]] constexpr size_t Dynamic = 0;

    enum class Backend {
        Base,
        MKL
    };

    __host__ __device__ consteval inline static bool IsHost() {
    #ifdef __CUDA_ARCH__
        return false;
    #else
        return true;
    #endif
    }

    [[maybe_unused]] __host__ __device__ consteval inline static bool IsDevice() {
        return !IsHost();
    }

    [[maybe_unused]] consteval inline static bool IsMSVC() {
    #ifdef _MSC_VER
        return true;
    #else
        return false;
    #endif
    }

    [[maybe_unused]] consteval inline static bool HasHDF5() {
    #ifdef PHYSICA_HDF5
        return true;
    #else
        return false;
    #endif
    }

    [[maybe_unused]] consteval inline static bool HasMKL() {
    #ifdef PHYSICA_MKL
        return true;
    #else
        return false;
    #endif
    }

    [[maybe_unused]] consteval inline static bool HasMPI() {
    #ifdef PHYSICA_MPI
        return true;
    #else
        return false;
    #endif
    }

    [[maybe_unused]] consteval inline static bool HasCUDA() {
    #ifdef PHYSICA_CUDA
        return true;
    #else
        return false;
    #endif
    }
}
