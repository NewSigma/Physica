/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <system_error>
#include "../SystemError.h"
#include "Physica/Core/Exception/CUDA/cuDSS.cuh"

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
        [[nodiscard]] const char* name() const noexcept final { return "cuDSS"; }
        [[nodiscard]] std::string message(int code) const final;
    };
}

std::string Impl::message(int code) const {
    switch (cudssStatus_t(code)) {
    case CUDSS_STATUS_SUCCESS:
        return "No error";
    case CUDSS_STATUS_NOT_INITIALIZED:
        return "Not initialized";
    case CUDSS_STATUS_ALLOC_FAILED:
        return "Allocation failed";
    case CUDSS_STATUS_INVALID_VALUE:
        return "Invalid value";
    case CUDSS_STATUS_NOT_SUPPORTED:
        return "Not supported";
    case CUDSS_STATUS_EXECUTION_FAILED:
        return "Execution failed";
    case CUDSS_STATUS_INTERNAL_ERROR:
        return "Internal error";
    case CUDSS_STATUS_IR_FAILED:
        return "Iterative refinement failed";
    default:
        return "Unknown error";
    }
}

void Physica::check(cudssStatus_t err) {
    constinit static const Impl instance{};
    if (err != CUDSS_STATUS_SUCCESS) [[unlikely]]
        throw SystemError(int(err), instance);
}
