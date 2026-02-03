/*
 * Copyright 2025-2026 Weibo He.
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

#include <benchmark/benchmark.h>
#include <source_location>

namespace Physica {
    std::string makeBenchID(std::source_location loc = std::source_location::current()) noexcept;
    size_t makeVectorSize(int64_t level, size_t sizeElem) noexcept;
    size_t makeMatrixSize(int64_t level, size_t sizeElem) noexcept;
}

#ifdef PHYSICA_LLVMIR
    #define PHYSICA_BENCH(x) [[clang::noinline]] x
#else
    #define PHYSICA_BENCH(x) x
#endif
