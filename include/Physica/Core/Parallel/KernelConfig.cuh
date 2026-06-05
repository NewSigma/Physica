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
#pragma once

#include <format>
#include "Physica/Macro.h"

namespace Physica {
    class KernelConfig {
        using This = KernelConfig;
    public:
        dim3 blocks;
        dim3 threads;
    public:
        __host__ __device__ constexpr KernelConfig(dim3 blocks_, dim3 threads_) noexcept : blocks(blocks_), threads(threads_) {}
        constexpr KernelConfig(const This&) = default;
        constexpr KernelConfig(This&&) noexcept = default;
        constexpr ~KernelConfig() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Operations */
        PHYSICA_API void dump();
    };
}

namespace std {
    template<>
    struct PHYSICA_API formatter<Physica::KernelConfig, char> {
        constexpr static auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
        static auto format(const Physica::KernelConfig& obj, std::format_context& ctx) -> std::format_context::iterator;
    };
}
