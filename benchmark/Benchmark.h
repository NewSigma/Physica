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

#include "Physica/Config.h"
#include <benchmark/benchmark.h>
#include <array>

namespace Physica {
    constexpr std::array<std::size_t, 6> CacheSizes{
        HostDevAttr::LineSizeL1D / 2,
        HostDevAttr::LineSizeL1D,
        HostDevAttr::CacheSizeL1D * 3 / 4,
        HostDevAttr::CacheSizeL2 * 3 / 4,
        HostDevAttr::CacheSizeL3 * 3 / 4,
        HostDevAttr::CacheSizeL3  * 6 / 5
    };
}

#ifdef PHYSICA_LLVMIR
    #define PHYSICA_BENCH(x) [[clang::noinline]] x
#else
    #define PHYSICA_BENCH(x) x
#endif
