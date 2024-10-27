/*
 * Copyright 2021 Weibo He.
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

#include <exception>
#include <Physica/Macro.h>

namespace Physica::Core {
    class PHYSICA_API IOException : public std::exception {
        using This = IOException;

        const char* msg;
    public:
        IOException(const char* msg_) : msg(msg_) {}
        IOException(const This&) = default;
        IOException(This&&) noexcept = default;
        ~IOException() noexcept override = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Getters */
        const char* what() const noexcept override { return msg; }
    };
}