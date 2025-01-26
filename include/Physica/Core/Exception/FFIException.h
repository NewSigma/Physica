/*
 * Copyright 2024-2025 Weibo He.
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
#include <system_error>
#include "Physica/Macro.h"

namespace Physica {
    class PHYSICA_API FFIException : public std::system_error {
        class Impl final : public std::error_category {
        public:
            [[nodiscard]] const char* name() const noexcept override final { return "FFI"; }
            [[nodiscard]] std::string message(int code) const override final {
                //using enum ffi_status; // TODO: clang 16 ~ 18 ICE, refactor once nvcc supports clang 19
                switch (code) {
                case ffi_status::FFI_OK:
                    return "No error";
                case ffi_status::FFI_BAD_ABI:
                    return "Bad ABI";
                case ffi_status::FFI_BAD_TYPEDEF:
                    return "Bad typedef";
                default:
                    return "Unknown";
                }
            }
        };
    public:
        FFIException(ffi_status code) : std::system_error(code, Impl()) {}
    };

    inline void check(ffi_status err) {
        if (err != FFI_OK) [[unlikely]]
            throw FFIException(err);
    }
}
