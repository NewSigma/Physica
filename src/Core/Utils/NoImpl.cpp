/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Utils/NoImpl.h"
#include <cassert>
#include <print>

namespace Physica {
    __host__ __device__ void noImpl(const char* msg) noexcept {
        assert(msg != nullptr);
    #ifdef __CUDA_ARCH__
        printf("[Error]: Not implemented\n%s\n", msg);
        __trap();
    #else
        std::println("[Error]: Not implemented\n{}", msg);
        std::abort();
    #endif
    }

    __host__ __device__ void noImpl(std::source_location loc) noexcept {
        noImpl(std::format("File: {}:{}:{}\nFunc: {}", loc.file_name(), loc.line(), loc.column(), loc.function_name()).c_str());
    }
}
