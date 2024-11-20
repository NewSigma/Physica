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

#include <system_error>
#include "Physica/Macro.h"
#ifdef PHYSICA_MKL
    #include <mkl_vsl.h>
#endif

namespace Physica::Core {
    class PHYSICA_API VSLException : public std::system_error {
        using Base = std::system_error;

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
            [[nodiscard]] const char* name() const noexcept override final { return "Intel MKL VSL"; }
            [[nodiscard]] std::string message(int) const override final;
        };
    public:
        VSLException(int code) noexcept : std::system_error(code, Impl()) {}
        /* Getters */
        [[nodiscard]] int code() const noexcept { return Base::code().value(); }
    };
}

namespace Physica {
    inline void vslCheckImpl(int err) {
        if (err != 0) [[unlikely]]
            throw Physica::Core::VSLException(err);
    }
}

#ifdef PHYSICA_MKL
    #define vslCheck(err) vslCheckImpl(err)
#else
    #define vslCheck(err)
#endif
