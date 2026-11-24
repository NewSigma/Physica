/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <system_error>
#include "../SystemError.h"
#include "Physica/Core/Exception/CUDA/cuRAND.cuh"

using namespace Physica;

namespace {
    class Impl final : public std::error_category {
    public:
        Impl() = default;
        Impl(const Impl&) = delete;
        Impl(Impl&&) noexcept = delete;
        ~Impl() = default;
        /* Operators */
        Impl& operator=(const Impl&) = delete;
        Impl& operator=(Impl&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] const char* name() const noexcept final { return "cuRAND"; }
        [[nodiscard]] std::string message(int code) const final;
    };
}

std::string Impl::message(int code) const {
    switch (code) {
    case CURAND_STATUS_SUCCESS:
        return "No error";
    case CURAND_STATUS_VERSION_MISMATCH:
        return "Header file and linked library version do not match";
    case CURAND_STATUS_NOT_INITIALIZED:
        return "Generator not initialized";
    case CURAND_STATUS_ALLOCATION_FAILED:
        return "Memory allocation failed";
    case CURAND_STATUS_TYPE_ERROR:
        return "Generator is wrong type";
    case CURAND_STATUS_OUT_OF_RANGE:
        return "Argument out of range";
    case CURAND_STATUS_LENGTH_NOT_MULTIPLE:
        return "Length requested is not a multple of dimension";
    case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED:
        return "GPU does not have double precision required by MRG32k3a";
    case CURAND_STATUS_LAUNCH_FAILURE:
        return "Kernel launch failure";
    case CURAND_STATUS_PREEXISTING_FAILURE:
        return "Preexisting failure on library entry";
    case CURAND_STATUS_INITIALIZATION_FAILED:
        return "Initialization of CUDA failed";
    case CURAND_STATUS_ARCH_MISMATCH:
        return "Architecture mismatch, GPU does not support requested feature";
    case CURAND_STATUS_INTERNAL_ERROR:
        return "Internal library error";
    default:
        return "Unknown error";
    }
}

void Physica::check(curandStatus_t err) {
    constinit static const Impl instance{};
    if (err != CURAND_STATUS_SUCCESS) [[unlikely]]
        throw SystemError(int(err), instance);
}
