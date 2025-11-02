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
    struct curandGenerator;
    using curandGenerator_t = curandGenerator*; // Do not conflict with other pointers 
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

#ifdef PHYSICA_MIMALLOC
    #include <mimalloc-override.h>
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

    enum class Backend : char {
        Base,
        MKL
    };

    __host__ __device__ consteval static bool IsHost() noexcept {
    #ifdef __CUDA_ARCH__
        return false;
    #else
        return true;
    #endif
    }

    [[maybe_unused]] __host__ __device__ consteval static bool IsDevice() noexcept {
        return !IsHost();
    }

    [[maybe_unused]] consteval static bool IsDebug() noexcept {
    #ifndef NDEBUG
        return true;
    #else
        return false;
    #endif
    }

    [[maybe_unused]] consteval static bool IsMSVC() noexcept {
    #ifdef _MSC_VER
        return true;
    #else
        return false;
    #endif
    }

    [[maybe_unused]] consteval static bool IsIntelLLVM() noexcept {
    #ifdef __INTEL_LLVM_COMPILER
        return true;
    #else
        return false;
    #endif
    }

    [[maybe_unused]] consteval static bool HasHDF5() noexcept {
    #ifdef PHYSICA_HDF5
        return true;
    #else
        return false;
    #endif
    }

    [[maybe_unused]] consteval static bool HasMKL() noexcept {
    #ifdef PHYSICA_MKL
        return true;
    #else
        return false;
    #endif
    }

    [[maybe_unused]] consteval static bool HasMPI() noexcept {
    #ifdef PHYSICA_MPI
        return true;
    #else
        return false;
    #endif
    }

    [[maybe_unused]] consteval static bool HasCUDA() noexcept {
    #ifdef PHYSICA_CUDA
        return true;
    #else
        return false;
    #endif
    }
}
