/*
 * Copyright 2021-2024 Weibo He.
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

#include <stdexcept>
#include <Physica/Macro.h>

namespace Physica::Core {
    class PHYSICA_API NoImplException : public std::runtime_error {
    public:
        constexpr static const char* DefaultMsg = "[Error]: Not implemented";
    public:
        NoImplException(const char* msg = DefaultMsg) : std::runtime_error(msg) {}
    };

    __host__ __device__ [[noreturn]] inline void noImpl(const char* msg = NoImplException::DefaultMsg) {
    #ifdef __CUDA_ARCH__
        printf("%s\n", msg);
        __trap();
    #else
        throw NoImplException(msg);
    #endif
    }
}
