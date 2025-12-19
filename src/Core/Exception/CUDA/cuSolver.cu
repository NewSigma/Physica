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
}

cuSolverException::cuSolverException(cusolverStatus_t code) noexcept : std::system_error(code, Impl()) {}
