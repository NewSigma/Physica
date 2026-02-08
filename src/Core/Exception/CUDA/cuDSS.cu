/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Exception/CUDA/cuDSS.cuh"
#include "Physica/Core/Utils/Builtin.h"

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

cuDSSException::cuDSSException(cudssStatus_t code) noexcept : std::system_error(code, Impl()) {}

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
    default:
        return "Unknown error";
    }
}
