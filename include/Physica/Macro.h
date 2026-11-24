/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2019-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <climits>
#include <cstddef>
#include <type_traits>
#include "Physica/Config.h" // IWYU pragma: export
/**
 * \file Macro.h: Non-Proliferation of Macros
 */
#ifdef PHYSICA_CUDA
    #include <cuda_runtime.h>
#else
    #define __host__
    #define __device__
#endif

#ifdef __GNUC__
    #define PHYSICA_API [[gnu::visibility("default")]]
    #define PHYSICA_API_GNU __attribute__((visibility("default"))) // More general than PHYSICA_API
#else
    #define PHYSICA_API __declspec(dllexport)
    #define PHYSICA_API_GNU
    #define __restrict__ __restrict

    #include <BaseTsd.h>
    using ssize_t = SSIZE_T;
#endif

namespace Physica {
    #ifndef PHYSICA_MKL
        using MKL_INT64 = long long int;
    #endif

    template<class T>
    class Traits;

    template<class T> requires(std::is_const_v<T> || std::is_volatile_v<T> || std::is_reference_v<T>)
    class Traits<T> : public Traits<std::remove_cvref_t<T>> {};
    /**
     * Forward declaration for friend-class-based tests
     */
    class Test;

    [[maybe_unused]] constexpr static unsigned int PhysicaWordSize = sizeof(void*) * CHAR_BIT;
    [[maybe_unused]] constexpr size_t Dynamic = 0;

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

    [[maybe_unused]] consteval static bool HasMimalloc() noexcept {
    #ifdef PHYSICA_MIMALLOC
        return true;
    #else
        return false;
    #endif
    }
}
