/*
 * Copyright 2024 Weibo He.
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

#include <ffi.h>
#include <Physica/Macro.h>
#include <stdexcept>

namespace Physica::Python {
    class PHYSICA_API FFIException : public std::exception {
        ffi_status status;
    public:
        FFIException(ffi_status status_) : status(status_) {}
        ~FFIException() noexcept override = default;
        inline const char* what() const noexcept override;
    };

    inline const char* FFIException::what() const noexcept {
        constexpr static const char* msg[3]{"Bad typedef", "Bad ABI", "Bad arg type"};
        return msg[int(status) - 1];
    }
}
