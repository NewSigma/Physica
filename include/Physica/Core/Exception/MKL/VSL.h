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

namespace Physica {
    class PHYSICA_API VSLException : public std::system_error {
        using Base = std::system_error;
    public:
        VSLException(int code) noexcept;
        /* Getters */
        [[nodiscard]] int code() const noexcept { return Base::code().value(); }
    };

    inline void check_vsl_impl(int err) {
        if (err != 0) [[unlikely]]
            throw VSLException(err);
    }
}

#ifdef PHYSICA_MKL
    #define check_vsl(err) check_vsl_impl(err)
#else
    #define check_vsl(err)
#endif
