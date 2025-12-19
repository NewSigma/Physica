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
#include <cudss.h>
#include "Physica/Macro.h"

namespace Physica {
    class PHYSICA_API cuDSSException : public std::system_error {
    public:
        cuDSSException(cudssStatus_t code) noexcept;
    };

    inline void check(cudssStatus_t err) {
        if (err != CUDSS_STATUS_SUCCESS) [[unlikely]]
            throw cuDSSException(err);
    }
}
