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
#pragma once

#include <system_error>
#include <cusolverDn.h>
#include "Physica/Macro.h"

namespace Physica {
    class PHYSICA_API cuSolverException : public std::system_error {
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
    public:
        cuSolverException(cusolverStatus_t code) noexcept : std::system_error(code, Impl()) {}
    };
}

namespace Physica {
    inline void check(cusolverStatus_t err) {
        if (err != CUSOLVER_STATUS_SUCCESS) [[unlikely]]
            throw cuSolverException(err);
    }
}
