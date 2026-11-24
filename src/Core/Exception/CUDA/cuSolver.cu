/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <system_error>
#include "../SystemError.h"
#include "Physica/Core/Exception/CUDA/cuSolver.cuh"

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
        [[nodiscard]] const char* name() const noexcept override final { return "cuSolver"; }
        [[nodiscard]] std::string message(int code) const override final;
    };
}

std::string Impl::message(int code) const {
    switch (code) {
    case CUSOLVER_STATUS_SUCCESS:
        return "No error";
    case CUSOLVER_STATUS_NOT_INITIALIZED:
        return "Not initialized";
    case CUSOLVER_STATUS_ALLOC_FAILED:
        return "Alloc failed";
    case CUSOLVER_STATUS_INVALID_VALUE:
        return "Invalid value";
    case CUSOLVER_STATUS_ARCH_MISMATCH:
        return "Arch mismatch";
    case CUSOLVER_STATUS_MAPPING_ERROR:
        return "Mapping error";
    case CUSOLVER_STATUS_EXECUTION_FAILED:
        return "Execution failed";
    case CUSOLVER_STATUS_INTERNAL_ERROR:
        return "Internal error";
    case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
        return "Matrix type not support";
    case CUSOLVER_STATUS_NOT_SUPPORTED:
        return "Not support";
    case CUSOLVER_STATUS_ZERO_PIVOT:
        return "Zero pivot";
    case CUSOLVER_STATUS_INVALID_LICENSE:
        return "Invalid license";
    case CUSOLVER_STATUS_IRS_PARAMS_NOT_INITIALIZED:
        return "IRS: Param not initialized";
    case CUSOLVER_STATUS_IRS_PARAMS_INVALID:
        return "IRS: Invalid param";
    case CUSOLVER_STATUS_IRS_PARAMS_INVALID_PREC:
        return "IRS: Invalid prec";
    case CUSOLVER_STATUS_IRS_PARAMS_INVALID_REFINE:
        return "IRS: Invalid refine";
    case CUSOLVER_STATUS_IRS_PARAMS_INVALID_MAXITER:
        return "IRS: Invalid matrix iter";
    case CUSOLVER_STATUS_IRS_INTERNAL_ERROR:
        return "IRS: Internal error";
    case CUSOLVER_STATUS_IRS_NOT_SUPPORTED:
        return "IRS: Not support";
    case CUSOLVER_STATUS_IRS_OUT_OF_RANGE:
        return "IRS: Out of range";
    case CUSOLVER_STATUS_IRS_NRHS_NOT_SUPPORTED_FOR_REFINE_GMRES:
        return "IRS: NRHS not support for refine GMRES";
    case CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED:
        return "IRS: Info not initialized";
    case CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED:
        return "IRS: Info not destroyed";
    case CUSOLVER_STATUS_IRS_MATRIX_SINGULAR:
        return "IRS: Singular matrix";
    case CUSOLVER_STATUS_INVALID_WORKSPACE:
        return "CuSolver Exception";
    default:
        return "Unknown error";
    }
}

void Physica::check(cusolverStatus_t err) {
    constinit static const Impl instance{};
    if (err != CUSOLVER_STATUS_SUCCESS) [[unlikely]]
        throw SystemError(int(err), instance);
}
